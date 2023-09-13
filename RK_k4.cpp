#include <iostream>
#include<fstream>
#include <cmath>
using namespace std;

double time_sum_1 = 7200;
double time_sum_2 = 28800;
double out_step_1 = 5;
double out_step_2 = 5;
//max delta per step
double delta_set = 1e-4;
double time_now = 0;
double out_time = 0;
double step = 0.1;
double delta_max = 100;

double temp_b = 298.15;
double temp_e = 353.15;
double temp_q = 353.15;
double temp_p = temp_e - temp_b;
double temp_i = temp_p / time_sum_1;
double temp0 = temp_b;

double MW = 18e-3;
double MA = 29e-3;
double Mw = 18;
double Ma = 29;
double Rg = 8.314472;
double rho_s = 2057;
double rho_l = 1000;
double re = 8.8e-3;
double le = 61.8e-3;
double Ve = 3.14159 * re * re * le;
double V1 = 1;
double V2 = V1 * 2.94;

double eps = 0.4;
double As = 6 * eps / (15e-6);
double ee = 2.71828;
double alp = 1e-3;
double kcou = 6.5e-5;

double c_w_s0 = 80;
double c_w0 = 0.025;
double c_a0 = 40.316;

double p_w = c_w0 * Rg * temp_e;
double p_a = c_a0 * Rg * temp_e;
double p_t = p_w + p_a;
double p_out = 1e5;

double KA = 10.074;
double KB = 1657.46;
double KC = -46.13;
double pse = pow(10, KA - (KB / (KC + temp_e)));
double cse = pse / Rg / temp_e;
double RHe = 0.02;

double p_w_e = pse * RHe;
double p_a_e = 101325 - p_w_e;

double aa = 1.017;
double bb = 1.678;
double cc = 5714.3;
double dd = 1.018;

double yw, ya, Aaw, Awa, dw = 0;
double lmdw, lmda, miuw, miua = 0;
double lmd, miu = 0;
double Tc = 0;

double p_star, c_star, Dg, Ev, Ev_b, hm, im, RHb, Xi, Xi_ppm, Xb = 0;
double c_w_st, c_wt, c_at, tempt = 0;
double c_w_si, c_wi, c_ai, tempi = 0;
double kc_w_s1, kc_w_s2, kc_w_s3, kc_w_s4 = 0;
double kc_w1, kc_w2, kc_w3, kc_w4 = 0;
double kc_a1, kc_a2, kc_a3, kc_a4 = 0;
double ktemp1, ktemp2, ktemp3, ktemp4 = 0;

double delta_c_w_s, delta_c_w, delta_c_a, delta_temp, delta_1, delta_2 = 0;

double dcs = 0;

int Jpl = 1800;
int Jtime = 73200 + Jpl;
double Jloua = 1;
double Jback = 1e-2;
double time_h = 0;

double dc_w_s(double c_w_s, double c_w, double c_a, double temp) {
    double p_t = (c_w + c_a) * Rg * temp;
    double Dg = 2.17e-5 * (101325 / p_t) * pow((temp / 273.15), 1.88);
    double hm = Dg * kcou;
    double p_star = pow(10, KA - (KB / (KC + temp)));
    double c_star = p_star / Rg / temp;
    double RHb = c_w / c_star;
    double Xb = 0.0002981 * pow(RHb, 0.2871);
    double Xi = c_w_s * MW / rho_s;
    double Evb = -Rg * temp * log(RHb);
    double Ev = Evb * aa * pow(ee, -bb * pow((cc * (Xi - Xb)), dd));
    double ic = hm * As * (pow(ee, -Ev / Rg / temp) * c_star - c_w);
    double d_c_w_s = -ic;
    return d_c_w_s;
}

double dc_w(double c_w_s, double c_w, double c_a, double temp) {
    double p_t = (c_w + c_a) * Rg * temp;
    double Dg = 2.17e-5 * (101325 / p_t) * pow((temp / 273.15), 1.88);
    double hm = Dg * kcou;
    double p_star = pow(10, KA - (KB / (KC + temp)));
    double c_star = p_star / Rg / temp;
    double RHb = c_w / c_star;
    double Xb = 0.0002981 * pow(RHb, 0.2871);
    double Xi = c_w_s * MW / rho_s;
    double Evb = -Rg * temp * log(RHb);
    double Ev = Evb * aa * pow(ee, -bb * pow((cc * (Xi - Xb)), dd));
    double ic = hm * As * (pow(ee, -Ev / Rg / temp) * c_star - c_w);
    double pmp = alp * (p_t - p_out);
    double d_c_w = (ic * V1 - alp * (p_t - p_out) * c_w) / (V1 + V2);
    return d_c_w;
}

double dc_a(double c_w, double c_a, double temp) {
    double p_t = (c_w + c_a) * Rg * temp;
    double pmp = alp * (p_t - p_out);
    double d_c_a = (-alp * (p_t - p_out) * c_a) / (V1 + V2);
    return d_c_a;
}

double dtemp(double temp) {
    double Q = (temp_q - temp) * lmd * 1 / 0.2e-3;
    double d_temp = Q / (2057 * 1700 * le * 1);
    return d_temp;
}

int main() {
    ofstream outfile1, outfile2;
    outfile1.open("F:/CppProject/RK_k4/out_1.txt");
    cout << p_w << "  " << p_a << endl;

    while (time_now <= time_sum_1) {
        step = 0.1;
        delta_max = 100;

        if (temp0 >= temp_e) {
            temp_q = temp_e;
        }

        while (delta_max > delta_set) {
            step = step / 2;

            Tc = temp0 - 273.15;
            yw = c_w0 / (c_w0 + c_a0);
            ya = 1 - yw;
            dw = c_w0 / c_a0;
            miuw = (-8e-7 * Tc * Tc + 4.01e-2 * Tc + 8.02) * 1e-6;
            miua = (-4e-6 * Tc * Tc + 4.81e-2 * Tc + 17.2) * 1e-6;
            lmdw = (1e-5 * Tc * Tc + 5.2e-3 * Tc + 1.83) * 1e-2;
            lmda = (-3e-7 * Tc * Tc + 7.7e-3 * Tc + 2.44) * 1e-2;
            Aaw = pow((pow(miua / miuw, 0.5) * pow(Mw / Ma, 0.25) + 1), 2) / pow(8 * (Ma / Mw + 1), 0.5);
            Awa = (miuw * Ma * Aaw) / (miua * Mw);
            miu = (pow(Ma, -0.5) * miua + dw * pow(Mw, -0.5) * miuw) / (pow(Ma, -0.5) + dw * pow(Mw, -0.5));
            lmd = ((ya * lmda) / (ya + yw * Aaw)) + ((yw * lmdw) / (yw + ya * Awa));

            c_w_si = c_w_s0;
            kc_w_s1 = dc_w_s(c_w_si, c_w0, c_a0, temp0);
            c_w_si = c_w_s0 + kc_w_s1 * step / 2;
            kc_w_s2 = dc_w_s(c_w_si, c_w0, c_a0, temp0);
            c_w_si = c_w_s0 + kc_w_s2 * step / 2;
            kc_w_s3 = dc_w_s(c_w_si, c_w0, c_a0, temp0);
            c_w_si = c_w_s0 + kc_w_s3 * step;
            kc_w_s4 = dc_w_s(c_w_si, c_w0, c_a0, temp0);

            c_wi = c_w0;
            kc_w1 = dc_w(c_w_s0, c_wi, c_a0, temp0);
            c_wi = c_w0 + kc_w1 * step / 2;
            kc_w2 = dc_w(c_w_s0, c_wi, c_a0, temp0);
            c_wi = c_w0 + kc_w2 * step / 2;
            kc_w3 = dc_w(c_w_s0, c_wi, c_a0, temp0);
            c_wi = c_w0 + kc_w3 * step;
            kc_w4 = dc_w(c_w_s0, c_wi, c_a0, temp0);

            c_ai = c_a0;
            kc_a1 = dc_a(c_w0, c_ai, temp0);
            c_ai = c_a0 + kc_a1 * step / 2;
            kc_a2 = dc_a(c_w0, c_ai, temp0);
            c_ai = c_a0 + kc_a2 * step / 2;
            kc_a3 = dc_a(c_w0, c_ai, temp0);
            c_ai = c_a0 + kc_a3 * step;
            kc_a4 = dc_a(c_w0, c_ai, temp0);

            tempi = temp0;
            ktemp1 = dtemp(tempi);
            tempi = temp0 + ktemp1 * step / 2;
            ktemp2 = dtemp(tempi);
            tempi = temp0 + ktemp2 * step / 2;
            ktemp3 = dtemp(tempi);
            tempi = temp0 + ktemp3 * step;
            ktemp4 = dtemp(tempi);

            c_w_st = c_w_s0 + (kc_w_s1 + 2 * kc_w_s2 + 2 * kc_w_s3 + kc_w_s4) * step / 6;
            c_wt = c_w0 + (kc_w1 + 2 * kc_w2 + 2 * kc_w3 + kc_w4) * step / 6;
            c_at = c_a0 + (kc_a1 + 2 * kc_a2 + 2 * kc_a3 + kc_a4) * step / 6;
            tempt = temp0 + (ktemp1 + ktemp2 * 2 + ktemp3 * 2 + ktemp4) * step / 6;

            dcs = (c_w_s0 - c_w_st) / step;

            delta_c_w_s = abs((c_w_st - c_w_s0) / c_w_s0);
            delta_c_w = abs((c_wt - c_w0) / c_w0);
            delta_c_a = abs((c_at - c_a0) / c_a0);
            delta_temp = abs((tempt - temp0) / temp0);
            delta_1 = max(delta_c_w, delta_c_a);
            delta_2 = max(delta_c_w_s, delta_temp);
            delta_max = max(delta_1, delta_2);
        }

        c_w_s0 = c_w_st;
        c_w0 = c_wt;
        c_a0 = c_at;
        temp0 = tempt;

        time_now = time_now + step;

        if (time_now >= out_time) {
            p_w = c_w0 * Rg * temp0;
            p_a = c_a0 * Rg * temp0;
            p_t = p_w + p_a;
            Xi = c_w_s0 * MW / rho_s;
            Xi_ppm = Xi * 1e6;
            out_time += out_step_1;
            time_h = time_now / 3600;

            outfile1 << time_h << "  " << p_w << "  " << p_a << "  " << Xi_ppm << " " << temp0 << " " << dcs << endl;
        }
    }

    cout << "1_done" << endl;

    kcou = 1.5e-5;
    c_w0 /= p_out;
    c_a0 /= p_out;
    p_out = 4e3;
    c_w0 *= p_out;
    c_a0 *= p_out;

    while (time_now <= time_sum_2) {
        step = 0.1;
        delta_max = 100;
        Jloua = 1;

        while (delta_max > delta_set) {
            step = step / 2;

            c_w_si = c_w_s0;
            kc_w_s1 = dc_w_s(c_w_si, c_w0, c_a0, temp0);
            c_w_si = c_w_s0 + kc_w_s1 * step / 2;
            kc_w_s2 = dc_w_s(c_w_si, c_w0, c_a0, temp0);
            c_w_si = c_w_s0 + kc_w_s2 * step / 2;
            kc_w_s3 = dc_w_s(c_w_si, c_w0, c_a0, temp0);
            c_w_si = c_w_s0 + kc_w_s3 * step;
            kc_w_s4 = dc_w_s(c_w_si, c_w0, c_a0, temp0);

            c_wi = c_w0;
            kc_w1 = dc_w(c_w_s0, c_wi, c_a0, temp0);
            c_wi = c_w0 + kc_w1 * step / 2;
            kc_w2 = dc_w(c_w_s0, c_wi, c_a0, temp0);
            c_wi = c_w0 + kc_w2 * step / 2;
            kc_w3 = dc_w(c_w_s0, c_wi, c_a0, temp0);
            c_wi = c_w0 + kc_w3 * step;
            kc_w4 = dc_w(c_w_s0, c_wi, c_a0, temp0);

            c_ai = c_a0;
            kc_a1 = dc_a(c_w0, c_ai, temp0);
            c_ai = c_a0 + kc_a1 * step / 2;
            kc_a2 = dc_a(c_w0, c_ai, temp0);
            c_ai = c_a0 + kc_a2 * step / 2;
            kc_a3 = dc_a(c_w0, c_ai, temp0);
            c_ai = c_a0 + kc_a3 * step;
            kc_a4 = dc_a(c_w0, c_ai, temp0);

            c_w_st = c_w_s0 + (kc_w_s1 + 2 * kc_w_s2 + 2 * kc_w_s3 + kc_w_s4) * step / 6;
            c_wt = c_w0 + (kc_w1 + 2 * kc_w2 + 2 * kc_w3 + kc_w4) * step / 6;
            c_at = c_a0 + (kc_a1 + 2 * kc_a2 + 2 * kc_a3 + kc_a4) * step / 6;

            dcs = (c_w_s0 - c_w_st) / step;

            delta_c_w_s = abs((c_w_st - c_w_s0) / c_w_s0);
            delta_c_w = abs((c_wt - c_w0) / c_w0);
            delta_c_a = abs((c_at - c_a0) / c_a0);

            delta_1 = max(delta_c_w, delta_c_a);
            delta_max = max(delta_1, delta_c_w_s);

        }

        c_w_s0 = c_w_st;
        c_w0 = c_wt;
        c_a0 = c_at;

        time_now = time_now + step;

        if (time_now >= out_time) {
            if (time_now >= Jtime) {
                p_w = p_w_e / 101325 * 4e3;
                p_a = p_w_e / 101325 * 4e3;
                c_w0 = p_w / Rg / temp0;
                c_a0 = p_a / Rg / temp0;

                Jtime += Jpl;
            }
            p_w = c_w0 * Rg * temp0;
            p_a = c_a0 * Rg * temp0;
            Xi = c_w_s0 * MW / rho_s;
            Xi_ppm = Xi * 1e6;
            out_time += out_step_1;
            time_h = time_now / 3600;

            outfile1 << time_h << "  " << p_w << "  " << p_a << "  " << Xi_ppm << " " << temp0 << " " << dcs << endl;
        }
    }

    cout << "2_done" << endl;
    outfile1.close();
    return 0;
}

