import numpy as np

def needleman_wunsch_linear(seq1, seq2, match=1, mismatch=-1, gap_linear=-4):
    n, m = len(seq1), len(seq2)
    # Матрица счета
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)
    # Матрица для прослеживания пути (0 - диагональ, 1 - сверху (делеция), 2 - слева (инсерция))
    traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # Инициализация нулевой строки и столбца
    for i in range(1, n + 1):
        score_matrix[i][0] = score_matrix[i-1][0] + gap_linear
        traceback_matrix[i][0] = 1  # Сверху (делеция в seq2, т.е. гэп в seq1)
    for j in range(1, m + 1):
        score_matrix[0][j] = score_matrix[0][j-1] + gap_linear
        traceback_matrix[0][j] = 2  # Слева (инсерция в seq2)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Совпадение или несовпадение
            diag_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            # Сдвиг вверх (гэп в seq2 -> seq1[i-1] напротив гэпа)
            up_score = score_matrix[i-1][j] + gap_linear
            # Сдвиг влево (гэп в seq1)
            left_score = score_matrix[i][j-1] + gap_linear

            # Находим максимальный счет
            max_score = max(diag_score, up_score, left_score)
            score_matrix[i][j] = max_score

            # Запоминаем, откуда пришли (для обратного хода)
            if max_score == diag_score:
                traceback_matrix[i][j] = 0
            elif max_score == up_score:
                traceback_matrix[i][j] = 1
            else:
                traceback_matrix[i][j] = 2

    return score_matrix, traceback_matrix

def needleman_wunsch_affine(seq1, seq2, match=1, mismatch=-1, gap_open=-10, gap_extend=-1):
    n, m = len(seq1), len(seq2)

    # Инициализация трех матриц для разных состояний
    M = np.zeros((n + 1, m + 1), dtype=int)  # Состояние Match/Mismatch (выравнивание буква-буква)
    X = np.zeros((n + 1, m + 1), dtype=int)  # Состояние Insertion (гэп в seq1 -> seq2 длиннее)
    Y = np.zeros((n + 1, m + 1), dtype=int)  # Состояние Deletion (гэп в seq2 -> seq1 длиннее)

    # Матрицы для обратного хода
    trace_M = np.zeros((n + 1, m + 1), dtype=int)
    trace_X = np.zeros((n + 1, m + 1), dtype=int)
    trace_Y = np.zeros((n + 1, m + 1), dtype=int)

    # Инициализация: бесконечно малые значения для недостижимых состояний
    INF = -10**9
    M[0][0] = 0
    X[0][0] = INF
    Y[0][0] = INF

    # Заполнение первого столбца (гэпы в seq2, т.е. seq2 пустая)
    for i in range(1, n + 1):
        M[i][0] = INF
        Y[i][0] = gap_open + (i-1) * gap_extend # Открываем гэп и продолжаем
        X[i][0] = INF
        trace_Y[i][0] = 2 # Все время из Y (продолжение гэпа)
        # Для первого шага (i=1) открытие гэпа
        if i == 1:
            Y[i][0] = gap_open

    # Заполнение первой строки (гэпы в seq1)
    for j in range(1, m + 1):
        M[0][j] = INF
        X[0][j] = gap_open + (j-1) * gap_extend
        Y[0][j] = INF
        trace_X[0][j] = 1
        if j == 1:
            X[0][j] = gap_open

    # Заполнение матриц
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # 1. Обновляем M (выравнивание i-го и j-го символов)
            # Приходим из M, X или Y + счет за совпадение/несовпадение
            from_M = M[i-1][j-1]
            from_X = X[i-1][j-1]
            from_Y = Y[i-1][j-1]
            best_prev = max(from_M, from_X, from_Y)
            M[i][j] = best_prev + (match if seq1[i-1] == seq2[j-1] else mismatch)
            # Запоминаем, откуда пришли в M
            if best_prev == from_M:
                trace_M[i][j] = 0
            elif best_prev == from_X:
                trace_M[i][j] = 1
            else:
                trace_M[i][j] = 2

            # 2. Обновляем X (гэп в seq1 -> двигаемся по seq2)
            # Приходим из M (открываем гэп) или из X (продолжаем гэп)
            from_M_to_X = M[i][j-1] + gap_open
            from_X_to_X = X[i][j-1] + gap_extend
            if from_M_to_X >= from_X_to_X:
                X[i][j] = from_M_to_X
                trace_X[i][j] = 0 # Пришли из M (открыли гэп)
            else:
                X[i][j] = from_X_to_X
                trace_X[i][j] = 1 # Пришли из X (продолжили гэп)

            # 3. Обновляем Y (гэп в seq2 -> двигаемся по seq1)
            # Приходим из M (открываем гэп) или из Y (продолжаем гэп)
            from_M_to_Y = M[i-1][j] + gap_open
            from_Y_to_Y = Y[i-1][j] + gap_extend
            if from_M_to_Y >= from_Y_to_Y:
                Y[i][j] = from_M_to_Y
                trace_Y[i][j] = 0 # Пришли из M (открыли гэп)
            else:
                Y[i][j] = from_Y_to_Y
                trace_Y[i][j] = 2 # Пришли из Y (продолжили гэп)

    # Финальный счет - максимум из трех состояний в правом нижнем углу
    final_score = max(M[n][m], X[n][m], Y[n][m])

    return M, X, Y, trace_M, trace_X, trace_Y, final_score

def traceback_linear(seq1, seq2, traceback_matrix):
    """Обратный ход для линейного штрафа."""
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback_matrix[i][j] == 0:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and traceback_matrix[i][j] == 1:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))

def traceback_affine(seq1, seq2, M, X, Y, trace_M, trace_X, trace_Y, final_score):
    """Обратный ход для аффинного штрафа с использованием трех матриц."""
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = len(seq1), len(seq2)

    # Определяем, в каком состоянии мы закончили
    if final_score == M[i][j]:
        state = 'M'
    elif final_score == X[i][j]:
        state = 'X'
    else:
        state = 'Y'

    while i > 0 or j > 0:
        if state == 'M':
            # Было выравнивание буква-буква
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            # Смотрим, откуда пришли в M
            if trace_M[i][j] == 0: # Из M
                state = 'M'
            elif trace_M[i][j] == 1: # Из X
                state = 'X'
            else: # Из Y
                state = 'Y'
            i -= 1
            j -= 1
        elif state == 'X':
            # Был гэп в seq1 (двигались по seq2) -> всталяем гэп в seq1
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            # Смотрим, откуда пришли в X
            if trace_X[i][j] == 0: # Из M (открыли гэп)
                state = 'M'
            else: # Из X (продолжили гэп)
                state = 'X'
            j -= 1
        else: # state == 'Y'
            # Был гэп в seq2 (двигались по seq1) -> всталяем гэп в seq2
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            # Смотрим, откуда пришли в Y
            if trace_Y[i][j] == 0: # Из M (открыли гэп)
                state = 'M'
            else: # Из Y (продолжили гэп)
                state = 'Y'
            i -= 1

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))

seq_a = "ATGCAGCAGCAGCCA"
seq_b = "ATATAT"

print("Алгоритм Нидлмана-Вунша")
print("\n--- Модель 1: Линейный штраф (Gap = -4) ---")
score_mat_lin, trace_mat_lin = needleman_wunsch_linear(seq_a, seq_b)
align1_lin, align2_lin = traceback_linear(seq_a, seq_b, trace_mat_lin)

print("Матрица счета (Linear):")
print(score_mat_lin)
print(f"\nФинальный Score: {score_mat_lin[len(seq_a)][len(seq_b)]}")
print("Выравнивание:")
print(align1_lin)
print(align2_lin)


print("\n--- Модель 2: Аффинный штраф (Open = -10, Extend = -1) ---")
M, X, Y, trace_M, trace_X, trace_Y, final_score_aff = needleman_wunsch_affine(seq_a, seq_b)
align1_aff, align2_aff = traceback_affine(seq_a, seq_b, M, X, Y, trace_M, trace_X, trace_Y, final_score_aff)

print("Матрица M (Match/Mismatch):")
print(M)
print("\nФинальный Score:", final_score_aff)
print("Выравнивание:")
print(align1_aff)
print(align2_aff)

print("Аффинная модель лучше описывает биологическую реальность возникновения инделов. " \
"Логично, что возниконовение мутационного события дорого, а его продление уже дешевеет как последствие первого")