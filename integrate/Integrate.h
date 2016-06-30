//
//  Integrate.h
//  integrate
//
//  Created by Дмитрий on 15.06.16.
//  Copyright © 2016 Дмитрий. All rights reserved.
//

#ifndef Integrate_h
#define Integrate_h
/*!
 \file
 \brief Заголовочный файл с описанием функций
 \author Dmitry
 \version 1.0
 Данный файл содержит в себе определения основных
 функций, используемых в демонстрационной программе
 */
typedef double(*function)(double);
typedef double(*dfunction)(double, double);
/*!
 Интегрирование методом прямоугольников
 \param[in] f Интегрируемая функция
 \param[in] a Начало интервала
 \param[in] b Конец интервала
 \param[in] step_count Количество шагов
 */
double integral_pram(function f, double a, double b, unsigned step_count);
///Дополнение для бесконечных пределов
double integral_pram_inf(function f, double a, double h, double eps);
/*!
 Интегрирование методом трапеций
 \param[in] f Интегрируемая функция
 \param[in] a Начало интервала
 \param[in] b Конец интервала
 \param[in] step_count Количество шагов
 */
double integral_trap(function f, double a, double b, unsigned step_count);
///Дополнение для бесконечных пределов
double integral_trap_inf(function f, double a, double h, double eps);
/*!
 Интегрирование методом Симпсона
 \param[in] f Интегрируемая функция
 \param[in] a Начало интервала
 \param[in] b Конец интервала
 \param[in] step_count Количество шагов
 */
double integral_simp(function f, double a, double b, unsigned step_count);
///Дополнение для бесконечных пределов
double integral_simp_inf(function f, double a, double h, double eps);
/*!
 Интегрирование методом Монте-Карло
 \param[in] f Интегрируемая функция
 \param[in] a Начало интервала
 \param[in] b Конец интервала
 \param[in] step_count Количество шагов
 */
double integral_monte(function f, double a, double b, unsigned step_count);
/*!
 Решение задачи Коши методом Рунге-Кутты 4ого порядка
 \param[out] x Выходной параметр
 \param[out] y Значение функции
 \param[in] f Входная функция
 \param[in] x0 Начало интервала
 \param[in] x1 Конец интервала
 \param[in] y0 Значение функции в точке x0
 \param[in] h Шаг сетки
 */
void integral_runge4(dfunction f, double x0, double x1, double y0, double* x, double* y, double h);
/*!
 Решение задачи Коши методом Рунге-Кутты 5ого порядка
 \param[out] x Выходной параметр
 \param[out] y Значение функции
 \param[in] f Входная функция
 \param[in] x0 Начало интервала
 \param[in] x1 Конец интервала
 \param[in] y0 Значение функции в точке x0
 \param[in] h Шаг сетки
 */
void integral_runge5(dfunction f, double x0, double x1, double y0, double* x, double* y, double h);
/*!
 Решение задачи Коши методом Рунге-Кутты 7ого порядка
 \param[out] x Выходной параметр
 \param[out] y Значение функции
 \param[in] f Входная функция
 \param[in] x0 Начало интервала
 \param[in] x1 Конец интервала
 \param[in] y0 Значение функции в точке x0
 \param[in] h Шаг сетки
 */
void integral_runge78(dfunction f, double x0, double x1, double y0, double* x, double* y, double h);
/*!
 Решение задачи Коши методом Эйлера
 \param[out] x Выходной параметр
 \param[out] y Значение функции
 \param[in] f Входная функция
 \param[in] x0 Начало интервала
 \param[in] x1 Конец интервала
 \param[in] y0 Значение функции в точке x0
 \param[in] h Шаг сетки
 */
void integral_eiler(dfunction f, double x0, double x1, double y0, double* x, double* y, double h);
/*!
 Апроксимация таблично заданной функции методом наименьших квадратов
 \param[out] a_res Массив выходных коэффициентов
 \param[in] n Количество заданных точек
 \param[in] x Массив значений x[i]
 \param[in] y Массив значений y[i]
 \param[in] s Степень апроксимирующего полинома
 */
void mnk(int s, int n, double *x, double *y, double *a_res);
//Вспомогательная функция возведения числа t в степень k
int power1(int t, int k);
#endif /* Integrate_h */