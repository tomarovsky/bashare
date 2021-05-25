#!/bin/bash

set -e -u

timeover=90
without_blocked_users=true
without_sbinnologin_users=true

print_usage() {
        echo "get sorted CSV file with users and occupied space"
        echo "Usage: '-l' your last login time over (default 90)"
	echo "-b flag to save blocked users to an outfile"
	echo "-n flag to save users with /sbin/nologin status"
}

while getopts 'l:bn' flag; do
        case "${flag}" in
                l) timeover="${OPTARG}" ;;
		b) without_blocked_users=false ;;
		n) without_sbinnologin_users=false ;;
                *) print_usage
                        exit 1 ;;
        esac
done

date_now=$(date +%Y%m%d)
output=users_over${timeover}_${date_now}.csv

# script
# get sort CSV file with users, dates and occupied space
ssh head02 lastlog -b ${timeover} |\
awk '!/Never log/ {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$4,$5,$6,$7,$8,$9); system("ssh meta01 rbh-du -H -u " $1 "| cut -f1")}' OFS='\t' 2> /dev/null | tail -n +2 |\
sort -t$'\t' -k8,8rh > ${output}
sed -i '/^root/d' ${output} # to remove root user

# condition for clearing the output file from blocked users
# condition for clearing users with the "/sbin/nologin" status
if [[ $without_blocked_users || $without_sbinnologin_users ]]; then
	touch tmp;
	cat ${output} | while read line; do
		username=`echo ${line} | awk '{print $1}'` ;
		if ${without_blocked_users}; then
			if [[ -z $(sudo passwd --status ${username} | grep "Password locked") ]]; then
				echo ${line} >> tmp;
			fi
		fi
		if ${without_sbinnologin_users}; then
			if [[ -z $(getent passwd ${username} | grep "/sbin/nologin") ]]; then
				echo ${line} >> tmp;
			fi
		fi
	done
	awk '!seen[$0]++' tmp > ${output};
	rm tmp;
fi

