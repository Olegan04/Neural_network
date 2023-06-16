#include <iostream>
#include <fstream>
#include <ostream>
#include <math.h>
#include <algorithm>
#include <sstream>
#include <ctime>
using namespace std;
double*** net;
double*** dnet;
double** zn;
double** znosh;
bool* ob;
int kolsl;
double** obdan;
int yk = 0, sr = 0, random = 55;
void Comput(int dl[]);
double Sigm(double x);
double PrSigm(int x, int y);
double Tang(double x);
double PrTang(int x, int y);
double ReLu(double x);
int PrReLU(int x, int y);
double LeakyReLu(double x);
double PrLR(int x, int y);
double PrHod(int dl[], int x, int y);
double Fun(int x, int y);
void OshExSl(int x, double y);
void OshSkrSl(int dl[], int x);
void Ves(int dl[], double a, double n);
void Obych(int dl[], double n1, double n2, int x, int y);
int POD(int par);
int main()
{   
    unsigned int start_time = clock();
    setlocale(LC_ALL, "Russian");
    fstream fin("Wes1.txt");
    fstream fn("Sloi.txt");
    fstream f("D:\\Документы\\Творческий проект\\Obrab\\Obrab\\Otv1.txt");
    string st;
    int kol = 0;
    int razm; double maxp, minp;
    f >> razm >> maxp >> minp;
    getline(f, st);
    getline(f, st);
    stringstream stream(st);
    while (getline(stream, st, ' ') ){
        kol++;
    }
    fn >> kolsl;
    int* dl = new int [kolsl];
    dl[0] = kol-1; 
    //Указываем количество неронов в слоях
    for (int i = 1; i < kolsl; i++) {
        fn >> dl[i];
    }
    ob = new bool[razm]; obdan = new double* [razm];
    net = new double** [kolsl - 1]; dnet = new double** [kolsl - 1]; zn = new double* [kolsl]; znosh = new double* [kolsl - 1];
    //Ваделяем память для масивов
    for (int i = 0; i < kolsl; i++) {
        if (i == 0)
            zn[0] = new double[dl[i]];
        else {
            net[i - 1] = new double* [dl[i]];
            dnet[i - 1] = new double* [dl[i]];
            for (int j = 0; j < dl[i]; j++) {
                net[i - 1][j] = new double[dl[i - 1] + 1];
                dnet[i - 1][j] = new double[dl[i - 1] + 1];
            }
            zn[i] = new double[dl[i]];
            znosh[i-1] = new double[dl[i]];
        }
        
    }
    //Считываем обучающие данные
    for (int i = 0; i < razm; i++) {
        obdan[i] = new double[kol];
        for (int j = 0; j < kol; j++) {
            f >> obdan[i][j];
            //cout << obdan[i][j] << ' ';
        }
        ob[i] = false;
        //cout << '\n';
    }
    f.close();
    //Считывание из файла веса нейронов
    for (int i = 0; i < kolsl-1; i++) {
        for (int j = 0; j < dl[i+1]; j++)
            for (int q = 0; q < dl[i] + 1; q++) {
                if (!fin.eof())
                    fin >> net[i][j][q];
                else net[i][j][q] = (rand() % 10 - 5) / 100.0;
                //cout << net[i][j][q]<<'\n';
                dnet[i][j][q] = 0;
            }
    }
    fin.close();
    f.close();
    fn.close();
    Obych(dl, 0.1, 0.1, razm, kol);
    unsigned int end_time = clock();
    unsigned int search_time = end_time - start_time;
    cout << '\n' << search_time / 1000.0;
    //Освобождаем память   
    for (int i = 0; i < kolsl-1; i++) {
        for (int j = 0; j < dl[i + 1]; j++)
            delete[] dnet[i][j];
        delete[] znosh[i];
        delete[] dnet[i];
    }
    for (int i = 0; i < razm; i++) {
        delete[] obdan[i];
    }
    delete[]ob;
    delete[]obdan;
    delete[]dnet;
    delete[]znosh;
    ofstream fout("Wes.txt");
    //Записываем новые веса в файл
    for (int i = 0; i < kolsl - 1; i++) {
        for (int j = 0; j < dl[i + 1]; j++)
            for (int q = 0; q < dl[i] + 1; q++)
                fout << net[i][j][q] <<'\n';
    }
    fout.close();
    cout << "\nОбучение завершенно\n";
   /* while (true)
    {
        for (int i = 0; i < kol-1; i++)
            cin >> zn[0][i];
        Comput(dl);
        cout << zn[kolsl - 1][0] << '\n';
    }*/
    fstream finpr("D:\\Документы\\Творческий проект\\Obrab\\Obrab\\Test.txt");
    while (!finpr.eof())
    {
        for (int i = 0; i < kol - 1; i++)
            finpr >> zn[0][i];
        Comput(dl);
        double protv;
        finpr >> protv;
        double x = zn[kolsl - 1][0] * (maxp - minp) + minp;
        double y = protv * (maxp - minp) + minp;
        cout << x << '\t' << y << '\t' << (max(protv, zn[kolsl - 1][0]) - min(protv, zn[kolsl - 1][0])) / max(protv, zn[kolsl - 1][0]) * 100 << "% \n";
    }
    return 0;
}
void Comput(int dl[]) {
    for (int i = 1; i < kolsl; i++)
        for (int j = 0; j < dl[i]; j++) {
            zn[i][j] = PrHod(dl, i - 1, j);
            //cout << zn[i][j] << '\n';
        }
}
double PrHod(int dl[], int x, int y) {
    double st = 0;
    for (int i = 1; i <= dl[x]; i++)
        st += zn[x][i - 1] * net[x][y][i];
    st += net[x][y][0];
    //cout << st << '\n';
    //cout << x + 2 <<' ' << kolsl << '\n';
    double neron;
    if (x+2 < kolsl)
        neron = Tang(st); 
    else
        neron = Sigm(st);
    //cout << neron << '\n';
    return neron;
}
double Sigm(double x) {
    return 1.0 / (1.0 + exp(-x));
}
double Tang(double x) {
    return (exp(x) - exp(-x)) / (exp(x) + exp(-x));
}
double ReLU(double x) {
    return max(0.0, x);
}
double LeakyReLu(double x) {
    if (x < 0) return 0.01 * x;
    return x;
}
void OshExSl(int x, double y) {
    for (int i = 0; i < x; i++) {
        znosh[kolsl - 2][i] = PrSigm(kolsl - 1, i) * (zn[kolsl - 1][i] - y);
        //cout << znosh[kolsl - 2][i] <<'\n';
    }
    
}
void OshSkrSl(int dl[], int x) {
    for (int i = x-1; i >= 0; i--)
        for (int j = 0; j < dl[i+1]; j++) {
            double sum = 0;
            for (int q = 0; q < dl[i + 2]; q++) {
                sum += znosh[i + 1][q] * net[i + 1][q][j+1];
                //cout << znosh[i + 1][q] * net[i + 1][q][j + 1] << '\n';
            }
                
            znosh[i][j] = PrTang(i+1, j) * sum;
            
           //cout << znosh[i][j]<<'\n';
        }
}
double PrSigm(int x, int y) {
    //cout << zn[x][y];
    return zn[x][y] * (1 - zn[x][y]);
}
double PrTang(int x, int y) {
    return 4 / pow(exp(zn[x][y]) + exp(-zn[x][y]), 2);
}
int PrReLu(int x, int y) {
    if (zn[x][y] >= 0) return 1;
    return 0;
}
double PrLR(int x, int y) {
    if (zn[x][y] >= 0) return 1;
    return 0.01;
}
void Ves(int dl[], double a, double n) {
    for(int i = kolsl-2; i >= 0 ; i--)
        for(int j = 0; j < dl[i+1]; j++)
            for (int q = 0; q < dl[i]+1; q++)
                if (q == 0) {
                    double x = a * znosh[i][j];
                    net[i][j][q] -= x +n * dnet[i][j][q];
                    dnet[i][j][q] = x;
                    //cout << net[i][j][q] << ' ' << '\n';
                }         
                else {
                    double x = a * znosh[i][j] * zn[i][q - 1];
                    net[i][j][q] -= x+n * dnet[i][j][q];
                    //cout << net[i][j][q] << ' ' << '\n';
                    dnet[i][j][q] = x;
                }                  
}
void Obych(int dl[], double n1, double n2, int x, int y) {
    int kol = 0, aaa = 0;
    double E = 0, SrE = 200.0, mizn = 10000.0;
    do
    {
        int i = POD(x);
        for (int j = 0; j < y - 1; j++)
            zn[0][j] = obdan[i][j];
        Comput(dl);
        
        OshExSl(dl[kolsl - 1], obdan[i][y-1]);
        OshSkrSl(dl, kolsl - 2);
        Ves(dl, n1, n2);
        for(int j = 0; j < dl[kolsl-1]; j++)
            E += 0.5 * (zn[kolsl-1][j] - obdan[i][y-1]) * (zn[kolsl-1][j] - obdan[i][y-1]);
        kol++;
        if (kol >= x) {
            aaa++;
            SrE = E / kol;
            if (SrE < mizn) mizn = SrE;
            E = 0;
            kol = 0;
            cout << '\n' << SrE;
        }
    } while (SrE > 3e-05);
    cout << "\nМинимум: " << mizn;
}
int POD(int par) {
    int x;
    random = (random + 5) * 7 / 4;
    srand(random);
    while (yk == sr) {
        x = rand() % par;
        if (ob[x] == false) {
            ob[x] = true;
            yk++;
        }
    }
    sr = yk;
    if (sr == par) {
        for (int i = 0; i < par; i++)
            ob[i] = false;
        sr = yk = 0;
    }
    return x;
}    
//0.000025 30 20 20
//0.000031 30 
//0.000031 30 20
//0.000032 20 20  20