#!/bin/bash

# Создаём папку fastqs
mkdir -p fastqs

# Цикл от 1 до 10
for i in {1..10}
do
    # Имя файла
    filename="fastqs/sample_${i}.fastq"
    
    # Записываем строку в файл
    echo "This is sample number $i" > "$filename"
done

echo "10 files created in fastqs/"
