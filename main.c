#include <stdio.h>
#include <math.h>

/*** МАТЕМАТИЧЕСКАЯ МОДЕЛЬ SEIR-D РАСПРОСТРАНЕНИЯ КОРОНАВИРУСА COVID 19 В НОВОСИБИРСКОЙ ОБЛАСТИ ***/

// N - вся популяция, S - восприимчивое население, E - бессимптомно инфицированное население, I - количество выявленных случаев заражения, R - количество вылечившиеся, D - количество умерших
#define MU 0.0188     // коэф. смертности от COVID-19
#define BETA 0.999    // скорость выздоровления заражённых случаев
#define RO 0.952      // скорость восстановления выявленных случаев
#define ALPHA_E 0.999 // коэф. заражения между бессимптомно инфицированным и восприимчивым населением
#define ALPHA_I 0.999 // коэф. заражения между инфицированным и восприимчивым населением (социальный фактор)
#define K 0.042       // частота появления симптомов в открытых случаях
#define N0 2798170    // всё население Новосибирской области
#define E0 99         // начальное количество бессимптомно инфицированных
#define R0 24         // начальное количество вылечившихся
#define GAMMA 0       // скорость повторного заражения (0 - устойчивый иммунитет)
#define C 1           // ограничение на передвижения граждан (изначально - 1 + C_ISOL * (...), сокращена до 1, т.к. C_ISOL = 0)

/******************* СИСТЕМА НЕЛИНЕЙНЫХ ОБЫКНОВЕННЫХ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ *******************/
double dS_dt(double S, double E, double I, double R, double D);
double dE_dt(double S, double E, double I, double R, double D);
double dI_dt(double S, double E, double I, double R, double D);
double dR_dt(double S, double E, double I, double R, double D);
double dD_dt(double S, double E, double I, double R, double D);
/************************************************************************************************ */

void methodEulerCauchy(double a, double b, double h, double *S, double *E, double *I, double *R, double *D);

int main() {
    int a = 0, b = 90;
    double eps = 1e-2;
    double h = 1;
    double e = E0, i = 0, r = R0, d = 0, s = N0 - i - e - r - d;
    printf("\n Исходные данные модели SEIR-D\n N0 = %d (всё население)\n S0 = %d (восприимчивое население)\n E0 = %d (бессимптомно инфицированные)\n I0 = %d (выявленные случаи / инфицированные с симптомами)\n R0 = %d (вылечившиеся)\n D0 = %d (умершие)\n a = %d (день начала отсчёта), b = %d (день конца отсчёта)\n\n", N0, (int)floor(s), E0, (int)floor(i), R0, (int)floor(d), a, b);

    methodEulerCauchy(a, b, h, &s, &e, &i, &r, &d);

    double delta, d1 = d;
    int k = 1;
    
    printf(" ННахождение решения системы...\n");
    do {
        h = h / 2;
        e = E0, i = 0, r = R0, d = 0, s = N0 - i - e - r - d;
        methodEulerCauchy(a, b, h, &s, &e, &i, &r, &d);
        delta = fabs(d - d1);
        d1 = d;
        printf(" %d: delta = %lf, h = %lf\n", k, delta, h);
        k++;
    } while (delta > eps);

    printf("\n Прогноз модели SEIR-D на %d день\n S = %d (восприимчивого населения)\n E = %d (бессимптомно инфицированных)\n I = %d (выявленных случаев / инфицированных с симптомами)\n R = %d (вылечившихся)\n D = %d (умерших)\n", b, (int)floor(s), (int)floor(e), (int)floor(i), (int)floor(r), (int)floor(d));
    double n = s + e + i + r + d;
    if ((int)round(n) == N0) {
        printf(" N = %d = N0 = %d\n Следовательно, всё население было учтено!\n\n", (int)round(n), N0);
    } else {
        printf(" Ошибка!\n N = %d != N0 = %d\n\n", (int)round(n), N0);
    }

    return 0;
}

double dS_dt(double S, double E, double I, double R, double D) {
    double N = S + E + I + R + D;
    return -C * (ALPHA_I * S * I + ALPHA_E * S * E) / N + GAMMA * R;
}

double dE_dt(double S, double E, double I, double R, double D) {
    double N = S + E + I + R + D;
    return C * (ALPHA_I * S * I + ALPHA_E * S * E) / N - (K + RO) * E;
}

double dI_dt(double S, double E, double I, double R, double D) {
    return K * E - BETA * I - MU * I + 0 * (S + R + D);
}

double dR_dt(double S, double E, double I, double R, double D) {
    return BETA * I + RO * E - GAMMA * R + 0 * (S + D);
}

double dD_dt(double S, double E, double I, double R, double D) {
    return MU * I + 0 * (S + E + R + D);
}

void methodEulerCauchy(double a, double b, double h, double *S, double *E, double *I, double *R, double *D) {
    int n = (int)ceil((b - a) / h) + 1;
    double s = *S, e = *E, i = *I, r = *R, d = *D;
    double si, ei, ii, ri, di;
    double s1, e1, i1, r1, d1;

    for (int k = 0; k <= n; k++) {
        s1 = s + h * dS_dt(s, e, i, r, d);
        e1 = e + h * dE_dt(s, e, i, r, d);
        i1 = i + h * dI_dt(s, e, i, r, d);
        r1 = r + h * dR_dt(s, e, i, r, d);
        d1 = d + h * dD_dt(s, e, i, r, d);
    
        si = s + (h / 2) * (dS_dt(s, e, i, r, d) + dS_dt(s1, e1, i1, r1, d1));
        ei = e + (h / 2) * (dE_dt(s, e, i, r, d) + dE_dt(s1, e1, i1, r1, d1));
        ii = i + (h / 2) * (dI_dt(s, e, i, r, d) + dI_dt(s1, e1, i1, r1, d1));
        ri = r + (h / 2) * (dR_dt(s, e, i, r, d) + dR_dt(s1, e1, i1, r1, d1));
        di = d + (h / 2) * (dD_dt(s, e, i, r, d) + dD_dt(s1, e1, i1, r1, d1));
        
        s = s1;
        e = e1;
        i = i1;
        r = r1;
        d = d1;
    }

    *S = s;
    *E = e;
    *I = i;
    *R = r;
    *D = d;
}
