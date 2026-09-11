#!/bin/bash
# Parallel FinaleMe hg19 reference downloader.
# Uses curl with HTTP Range for parallel byte-range requests on the big files.
set -u
cd /tmp/FinaleMe
mkdir -p data
DATA_DIR="$(pwd)/data"
echo "[$(date '+%H:%M:%S')] Starting parallel hg19 reference download"
echo "[$(date '+%H:%M:%S')] Free disk: $(df -h /Users/hermes | tail -1 | awk '{print $4}')"

# 1. Small files first (these finish in seconds)
download_small() {
    local url="$1"
    local dest="$2"
    local name="$(basename "$dest")"
    if [ -s "$dest" ]; then
        echo "[SKIP] $name ($(stat -f%z "$dest") bytes)"
        return 0
    fi
    echo "[GO] $name"
    curl -fSL --retry 5 --retry-delay 3 --connect-timeout 30 --max-time 600 \
        -o "$dest" "$url" 2>/dev/null
    if [ -s "$dest" ]; then
        echo "[DONE] $name ($(stat -f%z "$dest") bytes)"
    else
        echo "[FAIL] $name"
        rm -f "$dest"
    fi
}

# 2. Big file: parallel byte-range download using HTTP Range
download_big_parallel() {
    local url="$1"
    local dest="$2"
    local name="$(basename "$dest")"
    local total_size=750000000  # ~715 MB; we'll get the real size from headers

    if [ -s "$dest" ]; then
        local sz=$(stat -f%z "$dest")
        if [ "$sz" -gt 700000000 ]; then
            echo "[SKIP] $name ($sz bytes — looks complete)"
            return 0
        fi
        echo "[RESUME] $name has $sz bytes; continuing from offset"
    fi

    # Get the real size
    local real_size=$(curl -sIL --connect-timeout 30 "$url" | grep -i "^content-length" | tail -1 | awk '{print $2}' | tr -d '\r')
    if [ -z "$real_size" ]; then
        echo "[FAIL] $name: could not get content-length"
        return 1
    fi
    echo "[GO] $name: $real_size bytes from $url"

    # Split into 4 chunks
    local n_chunks=4
    local chunk_size=$((real_size / n_chunks))
    local tmp_prefix="${dest}.part."

    for i in $(seq 0 $((n_chunks - 1))); do
        local start=$((i * chunk_size))
        local end=$((start + chunk_size - 1))
        if [ $i -eq $((n_chunks - 1)) ]; then
            end=$((real_size - 1))
        fi
        local part="${tmp_prefix}${i}"
        if [ -s "$part" ]; then
            local psz=$(stat -f%z "$part")
            local expected=$((end - start + 1))
            if [ "$psz" -eq "$expected" ]; then
                continue  # chunk already done
            fi
        fi
        (
            curl -fSL --retry 5 --retry-delay 3 --connect-timeout 30 --max-time 1800 \
                --range "$start-$end" -o "$part" "$url" 2>/dev/null
        ) &
    done
    wait

    # Concatenate chunks
    cat "${tmp_prefix}"0 "${tmp_prefix}"1 "${tmp_prefix}"2 "${tmp_prefix}"3 > "$dest"
    rm -f "${tmp_prefix}"*

    local final_size=$(stat -f%z "$dest" 2>/dev/null || echo 0)
    if [ "$final_size" -eq "$real_size" ]; then
        echo "[DONE] $name ($final_size bytes)"
    else
        echo "[FAIL] $name: got $final_size, expected $real_size"
        rm -f "$dest"
    fi
}

# Step 1: parallel small files
echo ""
echo "=== Small files (Zenodo, fast) ==="
(
    download_small "https://zenodo.org/records/19392525/files/CG_motif.hg19.common_chr.pos_only.bedgraph.gz?download=1" "$DATA_DIR/CG_motif.hg19.common_chr.pos_only.bedgraph.gz"
) &
(
    download_small "https://zenodo.org/records/19392525/files/CpG_index.hg19.bed.gz?download=1" "$DATA_DIR/CpG_index.hg19.bed.gz"
) &
(
    download_small "https://zenodo.org/records/19392525/files/CpG_index.hg19.bed.gz.csi?download=1" "$DATA_DIR/CpG_index.hg19.bed.gz.csi"
) &
(
    download_small "https://zenodo.org/records/19392525/files/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed?download=1" "$DATA_DIR/wgEncodeDukeMapabilityRegionsExcludable_wgEncodeDacMapabilityConsensusExcludable.hg19.bed"
) &
(
    download_small "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes" "$DATA_DIR/hg19.chrom.sizes"
) &
wait

# Step 2: parallel chunks for big files
echo ""
echo "=== Big files (UCSC, 4-way parallel range) ==="
download_big_parallel "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit" "$DATA_DIR/hg19.2bit" &
hg19_pid=$!
download_big_parallel "https://zenodo.org/records/19392525/files/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw?download=1" "$DATA_DIR/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw" &
bw_pid=$!
wait $hg19_pid
wait $bw_pid

echo ""
echo "[$(date '+%H:%M:%S')] Done. Final state of $DATA_DIR:"
du -sh "$DATA_DIR"/*
ls -la "$DATA_DIR"
