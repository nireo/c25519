#!/usr/bin/env python3

from __future__ import annotations

import argparse
import io
import re
import sys
from pathlib import Path


LOCAL_INCLUDE_RE = re.compile(r'^\s*#include\s+"([^"]+)"\s*$')
SYSTEM_INCLUDE_RE = re.compile(r"^\s*#include\s+<([^>]+)>\s*$")
HEADER_GUARD_IFNDEF_RE = re.compile(r"^\s*#ifndef\s+([A-Za-z_][A-Za-z0-9_]*)\s*$")
HEADER_GUARD_DEFINE_RE = re.compile(r"^\s*#define\s+([A-Za-z_][A-Za-z0-9_]*)\s*$")
PP_IF_RE = re.compile(r"^\s*#\s*(if|ifdef|ifndef)\b")
PP_ENDIF_RE = re.compile(r"^\s*#\s*endif\b")
EXTERN_C_IFDEF_RE = re.compile(r"^\s*#ifdef\s+__cplusplus\s*$")
EXTERN_C_OPEN_RE = re.compile(r'^\s*extern\s+"C"\s*\{\s*$')
EXTERN_C_CLOSE_RE = re.compile(r'^\s*}\s*(?://\s*extern\s+"C")?\s*$')


def next_nonempty_line(lines: list[str], start: int) -> int | None:
    for index in range(start, len(lines)):
        if lines[index].strip():
            return index
    return None


def prev_nonempty_line(lines: list[str], start: int) -> int | None:
    for index in range(start, -1, -1):
        if lines[index].strip():
            return index
    return None


def strip_header_guard(lines: list[str]) -> list[str]:
    first = next_nonempty_line(lines, 0)
    if first is None:
        return lines

    second = next_nonempty_line(lines, first + 1)
    if second is None:
        return lines

    ifndef_match = HEADER_GUARD_IFNDEF_RE.match(lines[first])
    define_match = HEADER_GUARD_DEFINE_RE.match(lines[second])
    if not ifndef_match or not define_match:
        return lines
    if ifndef_match.group(1) != define_match.group(1):
        return lines

    last = prev_nonempty_line(lines, len(lines) - 1)
    if last is None or last <= second:
        return lines
    if not PP_ENDIF_RE.match(lines[last]):
        return lines

    return (
        lines[:first]
        + lines[first + 1 : second]
        + lines[second + 1 : last]
        + lines[last + 1 :]
    )


def strip_extern_c_wrapper(lines: list[str]) -> list[str]:
    nonempty = [index for index, line in enumerate(lines) if line.strip()]
    open_triplet: tuple[int, int, int] | None = None
    close_triplet: tuple[int, int, int] | None = None

    for i in range(len(nonempty) - 2):
        a, b, c = nonempty[i], nonempty[i + 1], nonempty[i + 2]
        if (
            EXTERN_C_IFDEF_RE.match(lines[a])
            and EXTERN_C_OPEN_RE.match(lines[b])
            and PP_ENDIF_RE.match(lines[c])
        ):
            open_triplet = (a, b, c)
            break

    for i in range(len(nonempty) - 2):
        a, b, c = nonempty[i], nonempty[i + 1], nonempty[i + 2]
        if (
            EXTERN_C_IFDEF_RE.match(lines[a])
            and EXTERN_C_CLOSE_RE.match(lines[b])
            and PP_ENDIF_RE.match(lines[c])
        ):
            close_triplet = (a, b, c)

    if open_triplet is None or close_triplet is None:
        return lines
    if open_triplet[0] >= close_triplet[0]:
        return lines

    removed = set(open_triplet + close_triplet)
    return [line for index, line in enumerate(lines) if index not in removed]


def expand_header(
    path: Path,
    root: Path,
    seen: set[Path],
    seen_unconditional_system_includes: set[str],
    stack: list[Path],
    out: io.StringIO,
) -> None:
    path = path.resolve()
    if path in seen:
        return
    if path in stack:
        chain = " -> ".join(p.name for p in stack + [path])
        raise ValueError(f"cyclic local include: {chain}")

    seen.add(path)
    stack.append(path)
    try:
        lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
        lines = strip_header_guard(lines)
        lines = strip_extern_c_wrapper(lines)
        conditional_depth = 0

        for line in lines:
            match = LOCAL_INCLUDE_RE.match(line)
            if match:
                include_path = (path.parent / match.group(1)).resolve()
                if not include_path.exists():
                    raise FileNotFoundError(f"missing local header: {include_path}")

                try:
                    include_path.relative_to(root)
                except ValueError as exc:
                    raise ValueError(
                        f"refusing to inline header outside repo: {include_path}"
                    ) from exc

                expand_header(
                    include_path,
                    root,
                    seen,
                    seen_unconditional_system_includes,
                    stack,
                    out,
                )
                continue

            match = SYSTEM_INCLUDE_RE.match(line)
            if match and conditional_depth == 0:
                include_name = match.group(1)
                if include_name in seen_unconditional_system_includes:
                    continue
                seen_unconditional_system_includes.add(include_name)

            out.write(line)

            if PP_IF_RE.match(line):
                conditional_depth += 1
                continue
            if PP_ENDIF_RE.match(line):
                conditional_depth -= 1
                if conditional_depth < 0:
                    raise ValueError(f"unbalanced preprocessor directives in {path}")
        if conditional_depth != 0:
            raise ValueError(f"unbalanced preprocessor directives in {path}")
    finally:
        stack.pop()


def build_single_header(input_path: Path, root: Path) -> str:
    out = io.StringIO()
    out.write("/* generated by make_single_header.py. do not edit directly. */\n")
    out.write("#ifndef C25519_SINGLE_H\n")
    out.write("#define C25519_SINGLE_H\n\n")
    out.write("#ifdef __cplusplus\n")
    out.write('extern "C" {\n')
    out.write("#endif\n\n")
    expand_header(input_path, root, set(), set(), [], out)
    out.write("\n#ifdef __cplusplus\n")
    out.write('} // extern "C"\n')
    out.write("#endif\n")
    out.write("\n#endif /* C25519_SINGLE_H */\n")
    return out.getvalue()


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent

    parser = argparse.ArgumentParser(
        description="Generate a single-header file for c25519."
    )
    _ = parser.add_argument(
        "output",
        nargs="?",
        default=str(repo_root / "c25519_single.h"),
        help="output file path, or '-' for stdout",
    )
    _ = parser.add_argument(
        "--input",
        default=str(repo_root / "c25519.h"),
        help="top-level header to expand",
    )
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    if not input_path.exists():
        parser.error(f"input header does not exist: {input_path}")
    try:
        input_path.relative_to(repo_root)
    except ValueError as exc:
        raise SystemExit(f"input header must be inside {repo_root}") from exc

    content = build_single_header(input_path, repo_root)

    if args.output == "-":
        sys.stdout.write(content)
        return 0

    output_path = Path(args.output).resolve()
    if output_path == input_path:
        parser.error("output path must differ from input path")
    if not output_path.parent.exists():
        parser.error(f"output directory does not exist: {output_path.parent}")

    output_path.write_text(content, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
