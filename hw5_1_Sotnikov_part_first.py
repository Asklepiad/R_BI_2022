#!/usr/bin/env python
# coding: utf-8

"""Copy of homework_lecture_5.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/10UA-JMol2zcM7kyKYwKzDuPftw3_tTx2

Всем привет! Пришло время нашей первой домашней работы не на степике. Надеюсь, будет весело :)

Так как, `pandas` вы еще не проходили, то я вам немного помогу. Эту домашку можно делать как в питоне, так и в R, само задание будет написано в `Google Colaboratory`.
"""

# Pandas понадобится нам для чтения денных
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# В переменную data_path надо положить путь до данных
data_path = "/home/asklepiad/bioinf/R/files/BI-stat-course/lecture_5_hypothesis_testing/homework/data"
expression_data = pd.read_csv(f"{data_path}/homework_lecture_5_data.csv", index_col=0)
expression_data.head()

b_cells_expression_data = expression_data.query("Cell_type == 'B_cell'")
nk_cells_expression_data = expression_data.query("Cell_type == 'NK_cell'")

"""В этом задании мы будем работать с данными об экспрессии генов в друх клеточных типах: в B-клетках и в NK-клетках. Выясним, средняя экспрессия каких генов значимо различается в этих клеточных типах.

Я буду показывать примеры на одном гене, а на основе них вы сможете сделать похожие задачи для всех генов.
"""

example_gene = "TMCC1"

"""## Задание 1

**2 баллов**

Посмотрим на распределение экспрессий гена `TMCC1` в обоих клеточных типах.
"""

sns.histplot(b_cells_expression_data[example_gene], stat="density");

sns.histplot(nk_cells_expression_data[example_gene], stat="density");

"""Кажется, что они немного различаются. Для начала давайте попробуем еще раз проверить центральную предельную теорему.

**Задание:**

Напишите функцию, которая будет принимать на вход экспрессии гена, семплировать их них выборки размера `sample_size`, считать среднюю экспрессию и повторять это `n_samples` раз. Примените эту функцию к экспрессиям гена `TMCC1` в обоих клеточных типах, визуализируйте их. Отличаются ли средние экспрессии данного гена у этих клеточных типов?
 
Сигнатура функции:

```python
def demonstrate_clt(expressions, sample_size, n_samples):
    mean_expressions = []
    for i in range(n_samples):
        sample = np.random.choice(expressions[example_gene], size=sample_size)
        mean_expressions.append(np.mean(sample))
    return mean_expressions
```
"""

def demonstrate_clt(expressions, sample_size, n_samples):
    mean_expressions = []
    for i in range(n_samples):
        sample = np.random.choice(expressions[example_gene], size=sample_size)
        mean_expressions.append(np.mean(sample))
    return mean_expressions

sns.histplot(b_means, stat="density")

sns.histplot(nk_means, stat="density")

"""А теперь посчитайте 95% доверительные интервалы для обоих распределений (примем тот факт, что средние распределены нормально для обоих клеточных типов) и скажите, отличается ли средняя экспрессия данного гена между клеточными типами?"""

b_mu = round(np.mean(b_cells_expression_data[example_gene]), 3)
nk_mu = round(np.mean(nk_cells_expression_data[example_gene]), 3)
b_std = np.std(b_cells_expression_data[example_gene])
nk_std = np.std(nk_cells_expression_data[example_gene])
b_se = b_std / np.sqrt(n_samples)
nk_se = nk_std / np.sqrt(n_samples)

f"The mean expression of {example_gene} gene in B-cells is {b_mu} with confidence interval [{round(b_mu-1.96*b_se, 3)}, {round(b_mu+1.96*b_se, 3)}]"

f"The mean expression of {example_gene} gene in NK-cells is {nk_mu} with confidence interval [{round(nk_mu-1.96*nk_se, 3)}, {round(nk_mu+1.96*nk_se, 3)}]"

"""Так как доверительные интервалы распределений пересекаются, можно предположить, что их они статически значимо не отличаются друг от друга. Однако, чтобы это утверждать, необходимо использовать статистические тесты (например, t-тест Стьюдента для несвязанных выборок в случае, если обе распределены нормально, или тест Манна-Уитни, если распределение выборок отличается от нормального(судя по диаграммам, это наш случай)).

## Задание 2

**4 баллов**

Вспомнили центральную предельную теорему и то, как считать доверительные интервалы в простом случае, теперь давайте воспользуемся библиотечной реализацией для того, чтобы протестировать уже все гены.
"""

import scipy.stats as st



"""Посчитаем доверительные интервалы для нашего демонстрационного гена в обоих клеточных типах:"""

# B клетки
st.t.interval(alpha=0.95, # 95% доверительный интервал
              df=len(b_cells_expression_data[example_gene]) - 1, # число степеней свободы - 1
              loc=np.mean(b_cells_expression_data[example_gene]), # Среднее
              scale=st.sem(b_cells_expression_data[example_gene])) # Стандартная ошибка среднего

# NK клетки
st.t.interval(alpha=0.95, # 95% доверительный интервал
              df=len(nk_cells_expression_data[example_gene]) - 1, # число степеней свободы - 1
              loc=np.mean(nk_cells_expression_data[example_gene]), # Среднее
              scale=st.sem(nk_cells_expression_data[example_gene])) # Стандартная ошибка среднего



"""Напишите функцию для проверки того, что доверительные интервалы пересекаются. На лекции мы тестировали гипотезы для равенства среднего выборки заданному числу и проверяли, попало ли оно в границы этого интервала или нет, если оно оказывалось за ними, то мы говорили, что средние отличаются. Здесь же мы имеем дело с двумя выборками, поэтому будем проверять, пересекаются ли доверительные интервалы, и, если нет, то говорить о том, что средние в выборках отличаются.

```python
def check_intervals_intersect(first_ci, second_ci):   
    if first_ci[0] < second_ci[0]:
        return second_ci[0] <= first_ci[1] # True or False
    elif first_ci[0] > second_ci[0]:
        return first_ci[0] <= second_ci[1] # True or False
    elif first_ci[0] == second_ci[0]:
        return True
```
"""



"""Теперь для каждого гена посчитайте доверительные интервалы в обоих клеточных типах, и проверьте, пересекаются ли они? Результаты можно добавлять в список, например:

```python
ci_test_results = []
for gene in genes:
    # B клетки
    b_ci = st.t.interval(confidence=0.95,
    df=len(b_cells_expression_data[gene]) - 1,
    loc=np.mean(b_cells_expression_data[gene]),
    scale=st.sem(b_cells_expression_data[gene])) 
    # NK клетки
    nk_ci = st.t.interval(confidence=0.95,
    df=len(nk_cells_expression_data[gene]) - 1,
    loc=np.mean(nk_cells_expression_data[gene]),
    scale=st.sem(nk_cells_expression_data[gene]))
    ci_test_results.append(check_intervals_intersect(b_ci, nk_ci))
```
"""



"""Попытайтесь оформить это в виде функции, которая будет принимать на вход две таблицы с экспрессиями и выдавать для каждого гена, значимо ли отличается его средняя экспрессия между клеточными типами.

```python
def check_dge_with_ci(first_table, second_table):
    # dge - differential gene expression

    return ci_test_results
```
"""

def check_dge_with_ci(first_table, second_table):
    # dge - differential gene expression
    genes1 = np.array(list(first_table.columns), dtype='str')[:-1]
    genes2 = np.array(list(second_table.columns), dtype='str')[:-1]
    for gene in genes1:
        if gene in genes2:
            ft_ci = st.t.interval(confidence=0.95,
            df = len(first_table[gene]) - 1,
            loc = np.mean(first_table[gene]),
            scale = st.sem(first_table[gene])) 
            # NK клетки
            sc_ci = st.t.interval(confidence=0.95,
            df = len(second_table[gene]) - 1,
            loc = np.mean(second_table[gene]),
            scale = st.sem(second_table[gene]))
            ci_test_results.append(check_intervals_intersect(ft_ci, sc_ci))
            print("Expression of gene", gene, check_intervals_intersect(ft_ci, sc_ci))
    return ci_test_results

"""## Задание 3

**4 баллов**

Давайте теперь применим для той же задачи `z-критерий`. Будем считать, что в данном случае $\alpha$ = 0.05, и если полученное `p-value` будет меньше, то экспрессия генов значимо отличается.
"""

from statsmodels.stats.weightstats import ztest

# Наш излюбленный ген
ztest(
    b_cells_expression_data[example_gene],
    nk_cells_expression_data[example_gene]
)



"""Попытайтесь оформить это в виде функции, которая будет принимать на вход две таблицы с экспрессиями и выдавать для каждого гена, значимо ли отличается его средняя экспрессия между клеточными типами.

```python
def check_dge_with_ztest(first_table, second_table):
    # dge - differential gene expression
    z_test_results = []
    genes1 = np.array(list(first_table.columns), dtype='str')[:-1]
    genes2 = np.array(list(second_table.columns), dtype='str')[:-1]
    for gene in genes1:
        if gene in genes2:
            ztest(first_table[gene], second_table[gene])
            z_test_results.append(ztest(first_table[gene], second_table[gene]))
    return z_test_results
```
"""



"""## Задание 4

**10 баллов**

Теперь пришла пора оформить все ваши старания в виде программы. Напишите программу, которая принимает на вход следующие аргуметры:

1. `first_cell_type_expressions_path` &ndash; путь до таблицы с экспрессиями генов для одного клеточного типа;
2. `second_cell_type_expressions_path` &ndash; путь до таблицы с экспрессиями генов для второго клеточного типа;
3. `save_results_table` &ndash; название таблицы с результатами.

Считывать аргументы можно любым удобным способом (например, `input`, `argparse`).

Как читать данные при помощи пандаса мы уже знаем, осталось понять, как записывать результаты. Допустим, вы записывали результаты ваших тестов в списки, тогда создать пандасовский датафрейм можно следующим образом:
"""

ci_test_results = [True, False, True]
z_test_results = [True, True, True]
# Опционально можно также сохранять p-value для z-критерия
z_test_p_values = [0.004, 0.01, 0.0001]
# Также сохраните разницу в средних экспрессиях между 1 и 2 таблицами для каждого гена,
# чтобы было понять, уменьшается или увеличивается экспрессия гена
mean_diff = [-10, 10, 0.5]

# Созданим словарь {'название колонки': список_значений}
results = {
    "ci_test_results": ci_test_results,
    "z_test_results": z_test_results,
    "z_test_p_values": z_test_p_values,
    "mean_diff": mean_diff
}

# Из словаря делаем датафрейм
results = pd.DataFrame(results)
results.head()

# Сохраним таблицу в .csv файл
results.to_csv("path_to_your_awesome_results.csv")

# Difexpression tool

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest


# Checking of confidence intervals intersecting
def check_intervals_intersect(first_ci, second_ci):   
    if first_ci[0] < second_ci[0]:
        return second_ci[0] <= first_ci[1] # True or False
    elif first_ci[0] > second_ci[0]:
        return first_ci[0] <= second_ci[1] # True or False
    elif first_ci[0] == second_ci[0]:
        return True
    
# Computing and checking CIs for two tables with genes   
def check_dge_with_ci(first_table, second_table, gene, ci_test_results):
    # dge - differential gene expression
    ft_ci = st.t.interval(confidence=0.95,
    df = len(first_table[gene]) - 1,
    loc = np.mean(first_table[gene]),
    scale = st.sem(first_table[gene])) 
    # NK клетки
    sc_ci = st.t.interval(confidence=0.95,
    df = len(second_table[gene]) - 1,
    loc = np.mean(second_table[gene]),
    scale = st.sem(second_table[gene]))
    ci_test_results.append(check_intervals_intersect(ft_ci, sc_ci))


# Computing z-score
def check_dge_with_ztest(first_table, second_table, gene, z_test_stat, z_test_pvalue):
    # dge - differential gene expression
    # ztest(first_table[gene], second_table[gene])
    #z_test_results.append(ztest(first_table[gene], second_table[gene]))
    z_test_stat.append(np.ravel(ztest(first_table[gene], second_table[gene]))[1]<0.05)
    z_test_pvalue.append(np.round(np.ravel(ztest(first_table[gene], second_table[gene]))[1], 3))

# Computing difference between means of each cell-type

def diff_mean(first_table, second_table, gene, diff_means):
    first_mu = round(np.mean(first_table[gene]), 3)
    second_mu = round(np.mean(second_table[gene]), 3)
    diff_means.append(second_mu - first_mu)


def difexpression_tool(first_cell_type_expressions_path, second_cell_type_expressions_path, save_results_table):
    first_table = pd.read_csv(first_cell_type_expressions_path)
    second_table = pd.read_csv(second_cell_type_expressions_path)
    
    genes1 = np.array(list(first_table.columns), dtype='str')[:-1]
    genes2 = np.array(list(second_table.columns), dtype='str')[:-1]
    genes = []
    ci_test_results = []
    z_test_stat = []
    z_test_pvalue = []
    diff_means = []
    for gene in genes1:
        if gene in genes2:
            
# First column in the final table - gene's name
            genes.append(gene)
    
    
# Second column in the final table - is CIs intersect each other
            check_dge_with_ci(first_table, second_table, gene, ci_test_results)

    
# Third and fourth column in the final table - z-score and p.value
            check_dge_with_ztest(first_table, second_table, gene, z_test_stat, z_test_pvalue)
            #z_test_results = np.matrix(check_dge_with_ztest(first_table, second_table, gene))
            #z_test_stat = np.ravel(z_test_results[:,1])<0.05
            #z_test_pvalue = np.round(np.ravel(z_test_results[:,1]), 3)


# Fifth column in the final table
            diff_mean(first_table, second_table, gene, diff_means)

# Writing output
    results = {
        "gene": genes,
        "ci_test_results": ci_test_results,
        "z_test_results": z_test_stat,
        "z_test_p_values": z_test_pvalue,
        "mean_diff": diff_means
    }
    results = pd.DataFrame(results)
    results[1:].to_csv(save_results_table)

"""## Задание 5

**0.5 баллов (доп.)**

В онлайне сложно знакомиться, а особенно сейчас. Созвонитесь с кем-то из других студентов и прикрипите сюда скрин вашего созвона. Можно коротко описать, о чем вы говорили)
"""
