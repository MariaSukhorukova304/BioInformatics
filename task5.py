import matplotlib.pyplot as plt
import numpy as np
import os

path_human = r"D:\BioInformatics\Homework3\Homo sapiens mitochondrion.fasta"
path_neanderthal = r"D:\BioInformatics\Homework3\неандертальца.fasta"

# Параметры для DotPlot
WINDOW_SIZE = 10
STRINGENCY = 8
STEP = 5

# --- ФУНКЦИЯ ДЛЯ ЧТЕНИЯ FASTA ---
def read_fasta_simple(filepath):
    sequence = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('>'):
                sequence.append(line)
        
    full_sequence = ''.join(sequence).upper()
    print(f"Загружено: {os.path.basename(filepath)} (длина: {len(full_sequence)} п.н.)")
    return full_sequence

# Загружаем последовательности
print("Загрузка последовательностей...")
seq_human = read_fasta_simple(path_human)
seq_neanderthal = read_fasta_simple(path_neanderthal)

print(f"\nДлина генома человека: {len(seq_human)} п.н.")
print(f"Длина генома неандертальца: {len(seq_neanderthal)} п.н.")
print("\nСтроим DotPlot")

# Создаем пустой график
plt.figure(figsize=(12, 10))

# Будем собирать координаты для точек
x_coords = []
y_coords = []

total_iterations = (len(seq_neanderthal) - WINDOW_SIZE) // STEP
current_iteration = 0

# Проходим по последовательностям с шагом STEP
for i in range(0, len(seq_neanderthal) - WINDOW_SIZE, STEP):
    current_iteration += 1
    window_neanderthal = seq_neanderthal[i:i+WINDOW_SIZE]
    
    for j in range(0, len(seq_human) - WINDOW_SIZE, STEP):
        window_human = seq_human[j:j+WINDOW_SIZE]
        
        # Считаем совпадения
        matches = 0
        for k in range(WINDOW_SIZE):
            if window_neanderthal[k] == window_human[k]:
                matches += 1
        
        if matches >= STRINGENCY:
            # Сохраняем координаты точки
            x_coords.append(i + WINDOW_SIZE//2)
            y_coords.append(j + WINDOW_SIZE//2)
    
    # Показываем прогресс
    if current_iteration % 200 == 0:
        print(f"  Прогресс: {current_iteration}/{total_iterations} (найдено {len(x_coords)} точек)")

print(f"Построение завершено! Найдено {len(x_coords)} совпадающих участков.")

plt.figure(figsize=(12, 10))

plt.scatter(x_coords, y_coords, 
           c='black', 
           marker='.', 
           s=1,  # очень маленькие точки
           alpha=0.5,  # полупрозрачные
           rasterized=True)  # растеризация для экономии памяти

# Настройка осей
plt.xlabel("Неандерталец (позиция в п.н.)", fontsize=12)
plt.ylabel("Современный человек (позиция в п.н.)", fontsize=12)
plt.title(f"DotPlot митохондриальных геномов\n(Окно: {WINDOW_SIZE}, Совпадений: ≥{STRINGENCY}, Шаг: {STEP})", 
          fontsize=14)

# Устанавливаем границы осей
plt.xlim(0, len(seq_neanderthal))
plt.ylim(0, len(seq_human))


max_len = max(len(seq_neanderthal), len(seq_human))
plt.plot([0, max_len], [0, max_len], 'r--', linewidth=1, alpha=0.5, label='Идеальное совпадение')
plt.legend()


plt.savefig("mitochondrial_dotplot_highres.png", dpi=300, bbox_inches='tight')
print("  Сохранен PNG (300 dpi)")

plt.tight_layout()
plt.show()


print("""
Наблюдения на графике:
------------------------
Основная диагональ: четкая непрерывная линия - последовательности коллинеарны
Нет перпендикулярных линий = нет крупных инверсий
Нет разрывов диагонали = нет крупных делеций/инсерций

Вывод.
Митохондриальные геномы неандертальца и современного человека 
демонстрируют высокую степень синтении без крупных структурных 
перестроек. Наблюдаются только точечные различия, что подтверждает 
консервативность структуры мтДНК у гоминид и их близкое 
эволюционное родство.
""")

# Дополнительная статистика
print("\nДополнительная статистика:")
print(f"  Длина генома человека: {len(seq_human)} п.н.")
print(f"  Длина генома неандертальца: {len(seq_neanderthal)} п.н.")
print(f"  Найдено совпадающих участков: {len(x_coords)}")
print(f"  Примерная плотность совпадений: {len(x_coords) / (len(seq_human) * len(seq_neanderthal) / 1e6):.2f} точек на млн. п.н.")