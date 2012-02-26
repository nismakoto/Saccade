#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

// 刻み幅
#define H 0.001
// 初期時間
#define T_start 0.0
// 終了時間
#define T_finish 4.0
// 視線の総移動量
#define Delta_g 2.0
// 時定数
#define T1 0.15
#define T2 0.012
#define TN 25
// 初期視線移動量
#define G_init 0.0
// 眼球運動の速度
#define V_init 20.0
// 神経積分器の発火頻度
#define N_init 0.5
// バースト細胞の発火頻度
#define B_init 5.0
// リセッタブル神経積分器の発火頻度
#define S_init 1.0

// ルンゲ・クッタ法
void d1(int i, double **solution, double(*saccade)(double, double, double, double, double))
{
        solution[1][i] = saccade(solution[0][0], solution[0][1],  solution[0][2],  solution[0][3],  solution[0][4]);
}

void d2(int i, double **solution, double(*saccade)(double, double, double, double, double))
{
        solution[2][i] = saccade(solution[0][0]+(H/2.0)*solution[1][0], solution[0][1]+(H/2.0)*solution[1][1], solution[0][2]+(H/2.0)*solution[1][2], solution[0][3]+(H/2.0)*solution[1][3], solution[0][4]+(H/2.0)*solution[1][4]);
}

void d3(int i, double **solution, double(*saccade)(double, double, double, double, double))
{
        solution[3][i] = saccade(solution[0][0]+(H/2.0)*solution[2][0], solution[0][1]+(H/2.0)*solution[2][1], solution[0][2]+(H/2.0)*solution[2][2], solution[0][3]+(H/2.0)*solution[2][3], solution[0][4]+(H/2.0)*solution[2][4]);
}

void d4(int i, double **solution, double(*saccade)(double, double, double, double, double))
{
        solution[4][i] = saccade(solution[0][0]+H*solution[3][0], solution[0][1]+H*solution[3][1], solution[0][2]+H*solution[3][2], solution[0][3]+H*solution[3][3], solution[0][4]+H*solution[3][4]);
}

// 重み付け
void weighting(int i, double **solution)
{
        solution[0][i] = solution[0][i]+H*(solution[1][i]+2.0*solution[2][i]+2.0*solution[3][i]+solution[4][i])/6.0;
}

// 脳幹のサッカードに関与する細胞群の応答を近似する関数
double cells(double x)
{
        if(x < 0.0)
        {
                x = -800*(1.0-exp(-fabs(x)/6.0));
        }
        else
        {
                x = 800*(1.0-exp(-fabs(x)/6.0));
        }

        return (x);
}

// g,v,n,b,s各微分方程式
double fg(double g, double v, double n, double b, double s)
{
        return(v);
}

double fv(double g, double v, double n, double b, double s)
{
        return(-v*(1.0/T1+1.0/T2)+(-g+n+(T1+T2)*b)/(T1*T2));
}

double fn(double g, double v, double n, double b, double s)
{
        return(-n/TN+b);
}

double fb(double g, double v, double n, double b, double s)
{
        return((-b+cells(Delta_g-s))*s);
}

double fs(double g, double v, double n, double b, double s)
{
        return(b);
}

int main()
{
        int i, j;
        // ルンゲ・クッタ法と各微分法定式の関数ポインタ配列
        double (*saccade[5])(double, double, double, double, double) = {fg, fv, fn, fb, fs};
        void (*rungekutta[4])(int, double **solution, double(*saccade)(double, double, double, double, double)) = {d1, d2, d3, d4};

        // 領域の動的確保
        // 行を作る
        double **solution = new double*[5];
        // 列を作る
        for(i=0; i < 5; i++)
        {
                solution[i] = new double[5];
        }

        // 初期値の代入
        double t = T_start;
        solution[0][0] = G_init; 
        solution[0][1] = V_init;
        solution[0][2] = N_init;
        solution[0][3] = B_init; 
        solution[0][4] = S_init;

        // ファイルを開く
        ofstream out("saccade.dat");

        for(t=0.0; t <= T_finish; t+=H)
        {
                for(i=0; i < 4; i++)
                {
                        for(j=0; j < 5; j++)
                        {
                                // ルンゲ・クッタ法を実行.
                                // 関数ポインタと配列によって順次実行
                                rungekutta[i](j, solution, saccade[j]); 
                        }
                }
                for(j=0; j < 5; j++)
                {
                        // 重み付け
                        weighting(j, solution);
                }
                out << t << " " << solution[0][0] << endl;
        }

        // ファイルを閉じる
        out.close();

        //領域の解放
        // 列を解放
        for(i=0; i < 5; i++)
        {
                delete [] solution[i];
        }
        // 行を解放
        delete [] solution;
        return 0;
}
