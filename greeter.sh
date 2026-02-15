#!/bin/bash

# Проверка количества аргументов
if [ $# -lt 3 ]; then
    echo "Error: Not enough arguments provided"
    exit 1
fi

# Присваиваем аргументы переменным
FIRST_NAME=$1
LAST_NAME=$2
GROUP=$3

# Выводим приветствие
echo "Welcome, $FIRST_NAME $LAST_NAME from group $GROUP!"
