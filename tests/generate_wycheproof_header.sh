#!/usr/bin/env sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)

JSON_PATH=${1:-"$SCRIPT_DIR/ed25519_test.json"}
OUT_PATH=${2:-"$SCRIPT_DIR/ed25519_wycheproof_vectors.h"}

if ! command -v jq >/dev/null 2>&1; then
    echo "error: jq is required to generate wycheproof vectors" >&2
    exit 1
fi

if [ ! -f "$JSON_PATH" ]; then
    echo "error: json file not found: $JSON_PATH" >&2
    exit 1
fi

total=$(jq -r '.numberOfTests' "$JSON_PATH")
included=$(jq -r '[.testGroups[] as $g | $g.tests[] | select((.sig|length)==128)] | length' "$JSON_PATH")
skipped=$((total - included))

{
    cat <<'HDR'
#ifndef C25519_ED25519_WYCHEPROOF_VECTORS_H
#define C25519_ED25519_WYCHEPROOF_VECTORS_H

typedef struct wycheproof_case {
    uint32_t tc_id;
    const char* pk_hex;
    const char* msg_hex;
    const char* sig_hex;
    int expected_valid;
} wycheproof_case;

static const wycheproof_case wycheproof_cases[] = {
HDR

    jq -r '.testGroups[] as $g
        | $g.publicKey.pk as $pk
        | $g.tests[]
        | select((.sig|length)==128)
        | [(.tcId|tostring), $pk, .msg, .sig, (if .result=="valid" then "1" else "0" end)]
        | @tsv' "$JSON_PATH" \
        | awk -F '\t' '{ printf "    { %su, \"%s\", \"%s\", \"%s\", %s },\n", $1, $2, $3, $4, $5 }'

    cat <<EOF2
};

static const size_t wycheproof_case_count = sizeof(wycheproof_cases) / sizeof(wycheproof_cases[0]);
static const size_t wycheproof_total_count = ${total}u;
static const size_t wycheproof_skipped_count = ${skipped}u;

#endif /* C25519_ED25519_WYCHEPROOF_VECTORS_H */
EOF2
} > "$OUT_PATH"

echo "wrote $OUT_PATH (included $included/$total vectors; skipped $skipped unsupported signature-length cases)"
