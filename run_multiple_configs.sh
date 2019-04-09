#!/bin/bash

source ~/refractor/absco/env/bin/activate

num_processors=$(grep -c ^processor /proc/cpuinfo)
absco_sw_dir=$(dirname $0)
num_configs=$#

# Round up
group_size=$(awk "BEGIN {printf(\"%.0f\", $num_configs / $num_processors)}")

run_script=$(mktemp)

echo "Creating run script with groups of size: $group_size for $num_configs configs on $num_processors processors"
group_count=1
for config_fn in $*; do
    log_fn=$(echo $config_fn | sed 's/\.ini$/.log/g')
    if [[ "$config_fn" == "$log_fn" ]]; then
        echo "Error, config and log filename the same. Will not overwrite config"
        exit 1
    fi
    
    log_dir=$(grep ^intdir $config_fn | awk '{print $3}')

    echo -n "mkdir -p $log_dir ; "
    echo -n "$absco_sw_dir/run_LBLRTM_ABSCO.py -e2e -i $config_fn -y > $log_dir/$log_fn 2>&1"
    if [[ $group_count == $group_size ]]; then
        # Start a new line
        group_count=1
        echo ""
    else
        group_count=$(expr $group_count '+' 1)
        echo -n " ; "
    fi
done > $run_script

# Rewrite script to group each line and send the set of tasks to the background
sed -i -e 's/^/\( /' -e 's/$/ \) \&/' $run_script

num_jobs=$(cat $run_script | wc -l)
echo "Launching $num_jobs jobs"
sh $run_script
rm $run_script
