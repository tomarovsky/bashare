#!/bin/bash

useroverfile=''
date_now=$(date +%Y%m%d)
dryrun='false'

print_usage() {
	echo "send messages to users"
	echo "Usage: '-i' your CSV file from user_overtimer.sh"
	echo "-d flag for dryrun"
}

while getopts 'i:d' flag; do
	case "${flag}" in
		i) useroverfile="${OPTARG}" ;;
		d) dryrun='true' ;;
		*) print_usage
			exit 1 ;;
	esac
done

mymail="andrey.tomarovsky@gmail.com"
dryrun_output="user_messages_${date_now}.txt"

# script
cat ${useroverfile} | while read line; do
	username=`echo $line | awk '{print $1}'` ;
	usermail=`getent passwd | grep "^${username}:" | awk -F ':' '{print $5}'` ;
	lastdata=`echo $line | awk '{print $2,$3,$4,$5,$6,$7}'`;
	usedspace=`echo $line | awk '{print $NF}'` ;
	message="Приветствую, ${username}.\n
	Обратил внимание, что вы не логинились на вычислительный кластер 'Макарьич' более трех месяцев (${lastdata}), при этом сейчас у нас хранится ~${usedspace} ваших данных. Я хотел бы уточнить их актуальность и необходимость в сохранении вашего аккаунта. В случае отсутствия ответа в течение недели, данные автоматически будут удалены.\n
	С уважением,
	младший администратор кластера 'Макарьич'" ;
	if ${dryrun}; then
		echo -e "---- message to ${username}:${usermail} ----" ;
		echo -e ${message} | tee -a ${dryrun_output} ;
	else
		echo -e ${message} | mail -s "Аккаунт на кластере 'Макарьич'" "tomasananan@gmail.com,${mymail}" ;
		#echo -e ${message} | mail -s "Аккаунт на кластере 'Макарьич'" "${usermail},${mymail}" ;
	fi
done
