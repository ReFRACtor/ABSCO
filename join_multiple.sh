#!/bin/bash

in_dir=$1
out_dir=$2

function usage {
    print "$0 <input_directory> <output_directory>"
}

function spectral_range {
    h5dump -d Extent_Ranges -y -A 0  $1 | awk '/DATA {/,/}/ { out=$0; sub(/^.*DATA.*$/, "", out); sub(/.*}/, "", out); sub(/^[ ]+/, "", out) ; sub(/,/, "", out); print out;}' | xargs
}

if [ -z "$in_dir" ]; then
    echo "No input directory supplied"
    usage
    exit 1
fi

if [ -z "$out_dir" ]; then
    print "No output directory supplied"
    usage
    exit 1
fi

# Create output directory if it does not already exist
mkdir -p $out_dir

mol_names=$(find $in_dir -name "*.nc" -exec basename {} \; | sed 's/_.*$//g' | sort | uniq | xargs)

echo -e "Found the following molecule names:\n$mol_names"

for mol in $mol_names; do
    mol_files=( $(find $in_dir -name "${mol}_*.nc" | sort) )
    num_files=${#mol_files[@]}
    echo -e "\nFound $num_files files for $mol:"
    for fn in ${mol_files[@]}; do
        echo -e "\t$fn"
    done

    first_fn=${mol_files[0]}
    last_fn=${mol_files[-1]}

    if [ $num_files -gt 1 ]; then
        # Parse out the beginning and ending wn range to put into the output filename
        first_range=($(spectral_range $first_fn))
        last_range=($(spectral_range $last_fn))
        new_range=$(printf "%0.5d-%0.5d" ${first_range[0]} ${last_range[-1]})

        first_base=$(basename $first_fn)
        first_parts=(${first_base//_/ })

        postfix_parts=${first_parts[@]:2}
        new_postfix=${postfix_parts// /_}

        output_fn="${mol}_${new_range}_${new_postfix}"

        $(dirname $0)/join_tables.py ${mol_files[@]} -o $out_dir/$output_fn
    else
        ln -svf $(readlink -f $first_fn) $out_dir/$(basename $first_fn)
    fi
done
