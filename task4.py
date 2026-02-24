class SeedAndExtend:
    def __init__(self, reference, query, k=4, X=2, match=1, mismatch=-1, gap_open=-10, gap_extend=-1):
        """
        Инициализация алгоритма Seed-and-Extend
        
        Параметры:
        reference: референсная последовательность (D)
        query: искомая последовательность (Q)
        k: длина k-мера
        X: порог X-drop
        match: награда за совпадение
        mismatch: штраф за несовпадение
        gap_open: штраф за открытие гэпа
        gap_extend: штраф за продолжение гэпа
        """
        self.reference = reference
        self.query = query
        self.k = k
        self.X = X
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        
        # Индекс базы данных (референса)
        self.index = self.build_index()
        
    def build_index(self):
        """
        Фаза 0: Построение индекса базы данных
        Разбиваем референс на все k-меры и создаем словарь {k-mer: [позиции]}
        """
        index = {}
        for i in range(len(self.reference) - self.k + 1):
            kmer = self.reference[i:i+self.k]
            if kmer in index:
                index[kmer].append(i)
            else:
                index[kmer] = [i]
        return index
    
    def find_seeds(self):
        """
        Фаза 1: Seeding
        Разбиваем запрос на k-меры и ищем совпадения в индексе
        """
        seeds = []
        query_kmers = {}
        
        # Разбиваем запрос на k-меры
        for i in range(len(self.query) - self.k + 1):
            kmer = self.query[i:i+self.k]
            query_kmers[kmer] = i
            
            # Ищем в индексе
            if kmer in self.index:
                for ref_pos in self.index[kmer]:
                    seeds.append({
                        'kmer': kmer,
                        'query_pos': i,
                        'ref_pos': ref_pos,
                        'length': self.k
                    })
        
        print(f"Найденные k-меры в запросе: {query_kmers}")
        print(f"Индекс базы данных: {self.index}")
        return seeds
    
    def calculate_score(self, seq1, seq2):
        """
        Вспомогательная функция для подсчета счета выравнивания
        """
        score = 0
        i, j = 0, 0
        in_gap = False
        
        while i < len(seq1) and j < len(seq2):
            if seq1[i] == '-' or seq2[j] == '-':
                # Обработка гэпов
                if not in_gap:
                    score += self.gap_open
                    in_gap = True
                else:
                    score += self.gap_extend
                
                if seq1[i] == '-':
                    j += 1
                else:
                    i += 1
            else:
                # Совпадение или несовпадение
                in_gap = False
                if seq1[i] == seq2[j]:
                    score += self.match
                else:
                    score += self.mismatch
                i += 1
                j += 1
        
        return score
    
    def extend_left(self, seed):
        """
        Расширение влево от seed
        """
        query_pos = seed['query_pos']
        ref_pos = seed['ref_pos']
        k = seed['length']
        
        # Начальный счет от seed (идеальное совпадение)
        current_score = k * self.match
        max_score = current_score
        best_left_extension = 0  # ИНИЦИАЛИЗИРУЕМ ЗДЕСЬ, ДО ЦИКЛА
        
        print(f"\n=== Расширение влево для seed '{seed['kmer']}' ===")
        print(f"Начальный счет: {current_score}")
        
        left_seq_ref = []
        left_seq_query = []
        
        # Расширяемся влево, пока есть куда и не сработал X-drop
        steps = 0
        while query_pos > 0 and ref_pos > 0:
            steps += 1
            query_pos -= 1
            ref_pos -= 1
            
            # Получаем символы
            q_char = self.query[query_pos]
            r_char = self.reference[ref_pos]
            
            # Добавляем в начало списков
            left_seq_ref.insert(0, r_char)
            left_seq_query.insert(0, q_char)
            
            # Вычисляем изменение счета
            if q_char == r_char:
                delta = self.match
            else:
                delta = self.mismatch
            
            current_score += delta
            
            # Обновляем максимум
            if current_score > max_score:
                max_score = current_score
                best_left_extension = steps  # Запоминаем шаг с максимумом
            
            print(f"Шаг {steps}: {r_char} vs {q_char} -> delta={delta}, cur={current_score}, max={max_score}")
            
            # Проверка X-drop
            if max_score - current_score >= self.X:
                print(f"X-drop сработал! max - cur = {max_score - current_score} >= {self.X}")
                break
        
        # Обрезаем по лучшей позиции
        if best_left_extension < steps:
            # Оставляем только те шаги, которые дали максимальный счет
            left_seq_ref = left_seq_ref[-best_left_extension:] if best_left_extension > 0 else []
            left_seq_query = left_seq_query[-best_left_extension:] if best_left_extension > 0 else []
        
        print(f"Лучшее расширение влево: {best_left_extension} шагов, max_score={max_score}")
        
        return ''.join(left_seq_ref), ''.join(left_seq_query), max_score
    
    def extend_right(self, seed):
        """
        Расширение вправо от seed
        """
        query_pos = seed['query_pos'] + seed['length']
        ref_pos = seed['ref_pos'] + seed['length']
        k = seed['length']
        
        # Начальный счет от seed (идеальное совпадение)
        current_score = k * self.match
        max_score = current_score
        best_right_extension = 0  # ИНИЦИАЛИЗИРУЕМ ЗДЕСЬ, ДО ЦИКЛА
        
        print(f"\n=== Расширение вправо для seed '{seed['kmer']}' ===")
        print(f"Начальный счет: {current_score}")
        
        right_seq_ref = []
        right_seq_query = []
        
        # Расширяемся вправо
        steps = 0
        while query_pos < len(self.query) and ref_pos < len(self.reference):
            steps += 1
            
            # Получаем символы
            q_char = self.query[query_pos]
            r_char = self.reference[ref_pos]
            
            # Добавляем в конец списков
            right_seq_ref.append(r_char)
            right_seq_query.append(q_char)
            
            # Вычисляем изменение счета
            if q_char == r_char:
                delta = self.match
            else:
                delta = self.mismatch
            
            current_score += delta
            
            # Обновляем максимум
            if current_score > max_score:
                max_score = current_score
                best_right_extension = steps  # Запоминаем шаг с максимумом
            
            print(f"Шаг {steps}: {r_char} vs {q_char} -> delta={delta}, cur={current_score}, max={max_score}")
            
            # Проверка X-drop
            if max_score - current_score >= self.X:
                print(f"X-drop сработал! max - cur = {max_score - current_score} >= {self.X}")
                break
            
            query_pos += 1
            ref_pos += 1
        
        # Обрезаем по лучшей позиции
        if best_right_extension < steps:
            # Оставляем только те шаги, которые дали максимальный счет
            right_seq_ref = right_seq_ref[:best_right_extension]
            right_seq_query = right_seq_query[:best_right_extension]
        
        print(f"Лучшее расширение вправо: {best_right_extension} шагов, max_score={max_score}")
        
        return ''.join(right_seq_ref), ''.join(right_seq_query), max_score
    
    def run(self):
        """
        Запуск полного алгоритма
        """
        print("=" * 60)
        print("АЛГОРИТМ SEED-AND-EXTEND")
        print("=" * 60)
        print(f"Референс: {self.reference}")
        print(f"Запрос: {self.query}")
        print(f"Параметры: k={self.k}, X={self.X}")
        print(f"match={self.match}, mismatch={self.mismatch}")
        print(f"gap_open={self.gap_open}, gap_extend={self.gap_extend}")
        print("=" * 60)
        
        # Поиск семян
        seeds = self.find_seeds()
        
        if not seeds:
            print("Семена не найдены!")
            return None
        
        print(f"\nНайдено семян: {len(seeds)}")
        
        best_alignment = None
        best_score = float('-inf')
        
        # Для каждого семени делаем расширение
        for i, seed in enumerate(seeds):
            print(f"\n{'=' * 40}")
            print(f"Обработка семени {i+1}: {seed}")
            print(f"{'=' * 40}")
            
            # Расширение влево
            left_ref, left_query, left_max = self.extend_left(seed)
            
            # Расширение вправо
            right_ref, right_query, right_max = self.extend_right(seed)
            
            # Формируем полное выравнивание
            ref_aligned = left_ref + seed['kmer'] + right_ref
            query_aligned = left_query + seed['kmer'] + right_query
            
            # Полный счет (учитываем seed + расширения)
            full_score = self.calculate_score(ref_aligned, query_aligned)
            
            print(f"\nИТОГОВОЕ ВЫРАВНИВАНИЕ для семени {i+1}:")
            print(f"Ref:  {ref_aligned}")
            print(f"Query: {query_aligned}")
            print(f"Score: {full_score}")
            
            if full_score > best_score:
                best_score = full_score
                best_alignment = {
                    'seed': seed,
                    'left_ref': left_ref,
                    'left_query': left_query,
                    'right_ref': right_ref,
                    'right_query': right_query,
                    'ref_aligned': ref_aligned,
                    'query_aligned': query_aligned,
                    'score': full_score,
                    'left_max': left_max,
                    'right_max': right_max
                }
        
        return best_alignment


# Запуск алгоритма
if __name__ == "__main__":
    reference = "CTAGGATCCAGGCATACGA"
    query = "GGATCCATTCATTA"
    
    # Создаем экземпляр алгоритма
    aligner = SeedAndExtend(reference, query, k=4, X=2)
    
    # Запускаем
    result = aligner.run()
    
    # Выводим финальный результат
    print("\n" + "=" * 60)
    print("ЛУЧШЕЕ ВЫРАВНИВАНИЕ")
    print("=" * 60)
    
    if result:
        print(f"Лучший k-мер: '{result['seed']['kmer']}'")
        print(f"Позиция в запросе: {result['seed']['query_pos']}")
        print(f"Позиция в референсе: {result['seed']['ref_pos']}")
        print(f"\nИтоговый Smax (глобальный максимум): {result['score']}")
        print(f"Scur для левого расширения: {result['left_max']}")
        print(f"Scur для правого расширения: {result['right_max']}")
        print(f"\nРасширение влево:")
        print(f"  Референс: {result['left_ref']}")
        print(f"  Запрос:   {result['left_query']}")
        print(f"\nРасширение вправо:")
        print(f"  Референс: {result['right_ref']}")
        print(f"  Запрос:   {result['right_query']}")
        print(f"\nИТОГОВОЕ ВЫРАВНИВАНИЕ:")
        print(f"Ref:  {result['ref_aligned']}")
        print(f"Query: {result['query_aligned']}")
        
        # Визуализация совпадений
        visual = []
        for r, q in zip(result['ref_aligned'], result['query_aligned']):
            if r == q:
                visual.append('|')
            else:
                visual.append(' ')
        print(f"      {''.join(visual)}")
    else:
        print("Выравнивание не найдено!")