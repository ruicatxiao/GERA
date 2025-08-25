#!/usr/bin/env bash

set -eo pipefail

for b in autocycler canu.sh canu_trim.py flye.sh genome_size_raven.sh lja.sh metamdbg.sh metamdbg_filter.py miniasm.sh necat.sh nextdenovo.sh raven.sh redbean.sh; do
    command -v $b >/dev/null 2>&1 || { echo 2>&1 "ERROR: $b not found, was Autocycler installed correctly?"; exit 1; }
done

PROG_NAME=$(basename $0 ".sh")
OUT=$PWD
THREADS=$(getconf _NPROCESSORS_ONLN)
COUNT=4
KMER=51
SIZE='AUTO'
POSSIBLE_ASSEMBLERS=("canu" "flye" "lja" "metamdbg" "miniasm" "necat" "nextdenovo" "raven" "redbean")

function usage {
        echo 2>&1
        echo "Usage: $PROG_NAME <reads> [<reads> ...] [options]"
        echo "  <reads>            Long reads to assemble in fastq(.gz) format"
        echo "  -o, --out          Output directory (default: $OUT)"
        echo "  -t, --threads      #threads to use (default: $THREADS)"
        echo "  -c, --count        #subsampled read sets to output (default: $COUNT)"
        echo "  -k, --kmer         K-mer size for De Bruijn graph (default: $KMER)"
        echo "  -s, --size         Genome size (default: $SIZE)"
        echo "  -a, --assemblers   Assemblers to use (default: all available)"
        echo "                     Possible assemblers: ${POSSIBLE_ASSEMBLERS[@]}"
        echo "                     Note: this argument MUST BE WRAPPED in quotes"
        echo
        exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--out) OUT="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -c|--count) COUNT="$2"; shift 2 ;;
        -k|--kmer) KMER="$2"; shift 2 ;;
        -s|--size) SIZE="$2"; shift 2 ;;
        -a|--assemblers) CHOSEN_ASSEMBLERS=(); read -ra CHOSEN_ASSEMBLERS <<< "$2"; shift 2 ;;
        -h|--help) usage ;;
        -*) echo "ERROR: Unknown option '$1'" >&2; usage ;;
        *) break ;;
    esac
done

if [[ $# -eq 0 ]]; then
    echo "ERROR: At least one read file is required." >&2
    usage
fi

READ_FILES=()
shopt -s extglob
for f in "$@"; do
    [ ! -f "$f" ] && { echo "ERROR: File does not exist: $f" >&2; exit 1; }
    [[ "$f" == *.@(fastq.gz|fq.gz|fastq|fq) ]] && READ_FILES+=($f) || echo "WARNING: $f doesn't have fq|fastq(.gz) extension" >&2
done
shopt -u extglob
[ ${#READ_FILES[@]} -eq 0 ] && { echo "ERROR: No valid read files supplied" >&2; exit 1; }

if [ ! ${#CHOSEN_ASSEMBLERS[@]} -eq 0 ]; then
    ASSEMBLERS=()
    for assembler in "${CHOSEN_ASSEMBLERS[@]}"; do
        found=0
        for possible in "${POSSIBLE_ASSEMBLERS[@]}"; do
            if [[ "$assembler" == "$possible" ]]; then
                found=1
                ASSEMBLERS+=($assembler)
                break
            fi
        done
        if [[ $found -eq 0 ]]; then
            echo "Error: Invalid assembler specified: '$assembler'" >&2
            echo "Possible assemblers are: ${POSSIBLE_ASSEMBLERS[@]}" >&2
            exit 1
        fi
    done
else
    ASSEMBLERS=$POSSIBLE_ASSEMBLERS
fi

[[ "$COUNT" =~ ^[1-9]$ ]] || { echo "ERROR: --count must be between 1 and 9 (inclusive)" >&2; exit 1; }

echo
echo "Starting $PROG_NAME:    $(date)"
echo "  - Valid read files:   ${#READ_FILES[@]}"
echo "  - Output directory:   $OUT"
echo "  - Threads:            $THREADS"
echo "  - Subsampling sets:   $COUNT"
echo "  - Graph K-mer:        $KMER"
echo "  - Genome size:        $SIZE"
echo "  - Assembling with:    ${ASSEMBLERS[@]}"
echo "------------------------------------------------------------------"
echo
mkdir -vp $OUT
echo

TABLE=${OUT}/metrics.tsv
autocycler table > $TABLE

for reads in "${READ_FILES[@]}"; do
    sample=$(basename "$reads" | sed -E 's/\.(fastq\.gz|fq\.gz|fastq|fq)$//')
    echo "Assembling $sample"
    echo "-----------------------------"
    echo

    sample_out=${OUT}/${sample}
    assemblies=${sample_out}/assemblies
    subsampled_reads=${sample_out}/subsampled_reads
    autocycler_out=${sample_out}/autocycler_out

    if [ $SIZE == 'AUTO' ]; then
        echo "Getting genome size with genome_size_raven.sh"
        echo
        genome_size=$(genome_size_raven.sh $reads $THREADS)
    else
        genome_size=$SIZE
    fi

    autocycler subsample --reads $reads --out_dir $subsampled_reads --genome_size $genome_size --count $COUNT

    mkdir -vp $assemblies
    for assembler in "${ASSEMBLERS[@]}"; do
        for i in $(seq $COUNT); do
            echo
            echo "Assembling set $i with $assembler  -----------"
            echo
            ${assembler}.sh ${subsampled_reads}/sample_0${i}.fastq ${assemblies}/${assembler}_0${i} $THREADS $genome_size
        done
    done
    rm -rf $subsampled_reads

    autocycler compress -i $assemblies -a $autocycler_out -t $THREADS --kmer $KMER
    autocycler cluster -a $autocycler_out
    for c in ${autocycler_out}/clustering/qc_pass/cluster_*; do
        autocycler trim -c "$c"
        autocycler resolve -c "$c"
    done

    autocycler combine -a $autocycler_out -i ${autocycler_out}/clustering/qc_pass/cluster_*/5_final.gfa
    cp ${autocycler_out}/consensus_assembly.fasta ${OUT}/${sample}.fasta
    autocycler table -a $autocycler_out -n $sample >> $TABLE
    echo
    echo "Finished assembling $sample"
    echo "-----------------------------"
    echo
done

echo "Finished $PROG_NAME:    $(date)"
exit 0