package ru.startandroid.mymethods;

import android.app.Activity;
import android.support.annotation.MainThread;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.text.TextUtils;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.TextView;
import static java.lang.Math.*;
public class MainActivity extends Activity {

    TextView a;
    TextView b;
    TextView eps;
    EditText enter_a;
    EditText enter_b;
    EditText enter_eps;
    Button btnPassiveSearch;
    Button btnDihotomySearch;
    Button btnGoldenSection;
    Button btnFibonacciSearch;
    Button btnTangentMethod;
    Button btnNRMethod;
    Button btnSecantMethod;
    TextView xmin;
    TextView ymin;
    TextView X1;
    TextView X2;
    TextView L;
    TextView alpha;
    EditText enter_X1;
    EditText enter_X2;
    EditText enter_L;
    EditText enter_alpha;
    Button btnSuccessiveDisplacementMethod;
    Button btnGradMethod;
    Button btnGradConstStep;
    Button btnMNGS;
    Button btnDivergentSeriesMethod;
    TextView Xmin;
    TextView Ymin;
    Button btnAcceleratedGradMethod;
    Button btnRavineMethod;
    Button btnNR2Method;
    Button btnSR1Method;
    Button btnConjugateGradMethod;
    TextView Time;
    Button btnPenaltyMethod;

    //@Override
    public int F(int n)
    {
        int x1=1, x2=1;
        for(int i=3;i<=n;i++)
        {
            int t=x1+x2;
            x1=x2;
            x2=t;
        }
        return x2;
    }
    public double f(double x) { return x*x+2*(x*log10(x/E)-2); }
    public double df(double x) { return 2*(x+log10(x/E)+1/(log(10))); }
    public double d2f(double x) { return 2*(1+1/(x*log(10))); }

    public double PassiveSearch(double a1, double b1, double eps)
    {
        long n = Math.round(Math.ceil((b1-a1)/eps));
        eps = (b1-a1)/n;
        double Xmin=a1, Ymin = Xmin*Xmin + 2*(Xmin*log10(Xmin/E)-2);
        for(int i=0; i<=n; i++)
        {
            double x = a1+i*eps;
            //System.out.printf("f[i]= %f, ymin= %f%n", f[i], Ymin);
            if(f(x)<=Ymin)
            {
                Ymin = f(x);
                Xmin = x;
            }
            else
            {
                break;
            }
        }
        return Xmin;
    }
    public double DihotomySearch(double a, double b, double eps)
    {
        double delta=eps/5, c, d;
        while(b-a >= eps)
        {
            c=(a+b)/2-delta/2;
            d=(a+b)/2+delta/2;
            //System.out.printf("c= %f, d= %f%n", c, d);
            if(f(c)<=f(d))
            {
                b=d;
            }
            else
            {
                a=c;
            }
        }
        return (b+a)/2;
    }
    public double GoldenSection(double a, double b, double eps)
    {
        double c=a+(b-a)*(3-sqrt(5))/2, d=a+(b-a)*(sqrt(5)-1)/2;
        while(b-a >= eps)
        {
            if(f(c)<=f(d))
            {
                b=d;
                d=c;
                c=a+(b-a)*(3-sqrt(5))/2;
            }
            else
            {
                a=c;
                c=d;
                d=a+(b-a)*(sqrt(5)-1)/2;
            }
        }
        return (b+a)/2;
    }
    public double FibonacciSearch(double a, double b, double eps)
    {
        int n=0;
        do { n++; }
        while(F(n+2)<(b-a)/eps);
        double c, d;
        c=a+(b-a)*F(n)/F(n+2);
        d=a+(b-a)*F(n+1)/F(n+2);
        while(n>=2)
        {
            n--;
            if(f(c)<=f(d))
            {
                b=d;
                d=c;
                c=a+(b-a)*F(n)/F(n+2);
            }
            else
            {
                a=c;
                c=d;
                d=a+(b-a)*F(n+1)/F(n+2);
            }

        }
        return (b+a)/2;
    }
    public double TangentMethod(double a, double b, double eps)
    {
        double Xmin=0, c=(f(b)-f(a)+df(a)*a-df(b)*b)/(df(a)-df(b));
        while (abs(df(c))>eps & (b-a)>eps)
        {
            if(df(c)==0)
            {
                Xmin=c;
                break;
            }
            else if (df(c)>0) { b=c; }
            else { a=c; }
            c=(f(b)-f(a)+df(a)*a-df(b)*b)/(df(a)-df(b));
            Xmin=c;
        }
        return Xmin;
    }
    public double NewtonRaphsonMethod(double a,double eps)
    {
        double a1=a+2*eps;
        while(abs(df(a))>eps | abs(a1-a)>eps)
        {
            a1=a;
            a=a-df(a)/d2f(a);
        }
        return a;
    }
    public double SecantMethod(double a, double b,double eps)
    {
        double a1=a, a2=b, a3=a2-df(a2)*(a1-a2)/(df(a1)-df(a2));
        while(abs(df(a3))>eps | abs(a1-a2)>eps)
        {
            a1=a2;
            a2=a3;
            a3=a2-df(a2)*(a1-a2)/(df(a1)-df(a2));
        }
        return a3;
    }
    public double f(double[] x)
    {
        return 4*x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]+6*x[0]-x[1]-2;
    }
    public double[] df(double[] x)
    {
        double[] y={8*x[0]-2*x[1]+6, 2*x[1]-2*x[0]-1};
        return y;
    }
    public double dfaMNGS(double[] x, double alp)
    {
        return alp*2*(2*df(x)[0]*df(x)[1]+df(x)[1]*df(x)[1]+4*df(x)[0]*df(x)[0])-8*df(x)[0]*x[0]
                +2*x[0]*df(x)[1]+2*x[1]*df(x)[0]-2*x[1]*df(x)[1]-6*df(x)[0]+df(x)[1];
    }
    public double d2faMNGS(double[] x, double alp)
    {
        return 2*(2*df(x)[0]*df(x)[1]+df(x)[1]*df(x)[1]+4*df(x)[0]*df(x)[0]);
    }
    public double absF(double[] x)
    {
        double sum=0;
        for(int i=0;i<x.length;i++)
        {
            sum = sum + x[i]*x[i];
        }
        return sqrt(sum);
    }
    public double fa1(double[] x, double alp)
    {
        return 4*x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]+6*x[0]-x[1]-2+4*alp*alp+alp*(8*x[0]-2*x[1]+6);
    }
    public double fa2(double[] x, double alp)
    {
        return 4*x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]+6*x[0]-x[1]-2+alp*alp+alp*(2*x[1]-2*x[0]-1);
    }
    public double dfa1(double[] x, double alp) { return 8*alp+8*x[0]-2*x[1]+6; }
    public double dfa2(double[] x, double alp) { return 2*alp+2*x[1]-2*x[0]-1; }
    public double d2fa1(double[] x, double alp) { return 8; }
    public double d2fa2(double[] x, double alp) { return 2; }
    public double[] SuccessiveDisplacementMethod(double[] Xo, double eps)
    {
        double a=-5, b=5;
        // методом золотого сечения
       double[] Xc={0,0},Xd={0,0},x={Xo[0]+10*eps,Xo[1]+10*eps};
        while(abs(f(x) - f(Xo))>eps) {
            x[0]=Xo[0];
            x[1]=Xo[1];
            for (int i = 0; i < 2; i++) {
                a=-10;
                b=10;
                double c = a + (b - a) * (3 - sqrt(5)) / 2, d = a + (b - a) * (sqrt(5) - 1) / 2;
                while (b - a >= eps*eps) {
                    Xc[0] = Xo[0];
                    Xc[1] = Xo[1];
                    Xd[0] = Xo[0];
                    Xd[1] = Xo[1];
                    Xc[i] = Xo[i] + c;
                    Xd[i] = Xo[i] + d;
                    if (f(Xc) <= f(Xd)) {
                        b = d;
                        d = c;
                        c = a + (b - a) * (3 - sqrt(5)) / 2;
                    } else {
                        a = c;
                        c = d;
                        d = a + (b - a) * (sqrt(5) - 1) / 2;
                    }
                }
                Xo[i] += (a + b) / 2;
            }
        }
        return Xo;
    }
    public double[] GradMethod(double[] Xo, double eps, double L)
    {
        double alp0=0.1,alp=alp0;
        double E=0.5;
        double[] x={0,0};
        while(absF(df(x))>eps) {
            for (int i = 0; i < x.length; i++) {
                x[i] = Xo[i] - alp * df(Xo)[i];
            }
            while (f(x) - f(Xo) > -alp * E * absF(df(Xo)) * absF(df(Xo))) {
                alp = alp * L;
                for (int i = 0; i < x.length; i++) {
                    x[i] = Xo[i] - alp * df(Xo)[i];
                }
            }
            for (int i = 0; i < x.length; i++)
            {
                x[i] = Xo[i] - alp * df(Xo)[i];
            }
            Xo=x;
            alp=alp0;
        }
        return x;
    }
    public double[] GradConstStep(double[] Xo, double eps, double alp)
    {
        double[] x=Xo;
        while(absF(df(x))>eps)
        {
            for (int i = 0; i < x.length; i++) {
                x[i] = Xo[i] - alp * df(Xo)[i];
            }
            Xo=x;
        }
        return x;
    }
    public  double[] MNGS(double[] Xo, double eps)
    {
        double alp = 0.5;
        double a1=alp+2*eps;
        double[] x={0,0};
        while(absF(df(x))>eps)
        {
            while (abs(dfaMNGS(Xo, alp)) > eps*eps | abs(a1 - alp) > eps*eps) {
                a1 = alp;
                alp = alp - dfaMNGS(Xo, alp) / d2faMNGS(Xo, alp);
            }
            x[0] = Xo[0] - alp * df(Xo)[0];
            x[1] = Xo[1] - alp * df(Xo)[1];
            Xo[0]=x[0];
            Xo[1]=x[1];
        }
        return x;
    }
    public double[] DivergentSeriesMethod(double[] Xo, double eps)
    {
        int k=1;
        double[] x=Xo;
        while(absF(df(x))>eps)
        {
            for (int i = 0; i < x.length; i++) {
                x[i] = Xo[i] - df(Xo)[i]/k;
            }
            k++;
            Xo=x;
        }
        return x;
    }
    public double Fxalp(double[] A, double[] B, double alp)
    {
        // x=A+alp*B; A,B - vectors

        return 4*(A[0]+alp*B[0])*(A[0]+alp*B[0])+(A[1]+alp*B[1])*(A[1]+alp*B[1])-2*(A[0]+alp*B[0])*(A[1]+alp*B[1])+
                6*(A[0]+alp*B[0])-(A[1]+alp*B[1])-2;
    }
    public double[] AcceleratedGradMethod(double[] Xo, double eps)
    {
        double[] x=Xo;
        //System.out.printf("METHOD");
        double alp, a1, a, b, c, d;
        double[] y={0,0}, B;
        boolean flag = true;
        while(absF(df(x))>eps) {
            //----------------------------mngs: 2 steps
            alp = 1 / 2;
            a1 = alp + 2 * eps;
            for (int i = 0; i < 2; i++) {
                while (abs(dfaMNGS(x, alp)) > eps*eps | abs(a1 - alp) > eps*eps) {
                    a1 = alp;
                    alp = alp - dfaMNGS(x, alp) / d2faMNGS(x, alp);
                }
                y[0] = x[0] - alp * df(x)[0];
                y[1] = x[1] - alp * df(x)[1];
                x[0] = y[0];
                x[1] = y[1];
            }

            //---------------------------f(x)->min(alp) by Golden Section Method
            flag = true;
            a = -10;
            b = 10;
            double a2=a,b2=b;
            B=x;
            while (flag) {
                c = a + (b - a) * (3 - sqrt(5)) / 2;
                d = a + (b - a) * (sqrt(5) - 1) / 2;
                while (b - a >= eps*eps) {
                    B[0]=y[0]-x[0];
                    B[1]=y[1]-x[1];
                    if (Fxalp(x, B, c) <= Fxalp(x, B, d)) {
                        b = d;
                        d = c;
                        c = a + (b - a) * (3 - sqrt(5)) / 2;
                    } else {
                        a = c;
                        c = d;
                        d = a + (b - a) * (sqrt(5) - 1) / 2;
                    }
                }
                if ((a + b) / 2 >= a2 & (a + b) / 2 <= a2+eps) {
                    b = a + 2 * eps;
                    a = a - b2+a2;
                    a2=a;
                    b2=b;
                } else if ((a + b) / 2 >= b2-eps & (a + b) / 2 <= b2) {
                    a = b - 2 * eps;
                    b = b + b2-a2;
                    a2=a;
                    b2=b;
                } else {
                    flag = false;
                }
            }
            alp = (a + b) / 2;
            x[0] = x[0] + alp * (y[0] - x[0]);
            x[1] = x[1] + alp * (y[1] - x[1]);
        }
        return x;
    }
    public double[] RavineMethod(double[] Xo, double eps)
    {
        double[] Xo1={Xo[0]+eps/2, Xo[1]+eps/2},y={0,0},y1={0,0},B={0,0};
        double alp, a1, a, b, c, d;
        boolean flag;
        while(absF(df(Xo))>eps) {
            // MNGS(Xo)
            alp = 0.5;
            a1 = alp + 2 * eps;
            y[0] = 0;
            y[1] = 0;
            while (abs(dfaMNGS(Xo, alp)) > eps | abs(a1 - alp) > eps) {
                a1 = alp;
                alp = alp - dfaMNGS(Xo, alp) / d2faMNGS(Xo, alp);
            }
            y[0] = Xo[0] - alp * df(Xo)[0];
            y[1] = Xo[1] - alp * df(Xo)[1];
            Xo[0] = y[0];
            Xo[1] = y[1];
            // MNGS(Xo1)
            alp = 0.5;
            a1 = alp + 2 * eps;
            y1[0] = 0;
            y1[1] = 0;
            while (abs(dfaMNGS(Xo1, alp)) > eps*eps | abs(a1 - alp) > eps*eps) {
                a1 = alp;
                alp = alp - dfaMNGS(Xo1, alp) / d2faMNGS(Xo1, alp);
            }
            y1[0] = Xo1[0] - alp * df(Xo1)[0];
            y1[1] = Xo1[1] - alp * df(Xo1)[1];
            Xo1[0] = y1[0];
            Xo1[1] = y1[1];
            //---------------------------f(x)->min(alp) by Golden Section Method
            flag = true;
            a = -10;
            b = 10;
            double a2=a,b2=b;
            while (flag) {
                c = a + (b - a) * (3 - sqrt(5)) / 2;
                d = a + (b - a) * (sqrt(5) - 1) / 2;
                while (b - a >= eps*eps) {
                    B[0] = y[0] - y1[0];
                    B[1] = y[1] - y1[1];
                    if (Fxalp(y, B, c) <= Fxalp(y, B, d)) {
                        b = d;
                        d = c;
                        c = a + (b - a) * (3 - sqrt(5)) / 2;
                    } else {
                        a = c;
                        c = d;
                        d = a + (b - a) * (sqrt(5) - 1) / 2;
                    }
                }
                if ((a + b) / 2 >= a2 & (a + b) / 2 <= a2+eps) {
                    b = a + 2 * eps;
                    a = a - b2+a2;
                    a2=a;
                    b2=b;
                } else if ((a + b) / 2 >= b2-eps & (a + b) / 2 <= b2) {
                    a = b - 2 * eps;
                    b = b + b2-a2;
                    a2=a;
                    b2=b;
                } else {
                    flag = false;
                }
            }
            alp = (a + b) / 2;
            Xo[0] = y[0] + alp * (y[0] - y1[0]);
            Xo[1] = y[1] + alp * (y[1] - y1[1]);
        }

        return Xo;
    }

    public double[] NR2Method(double[] Xo, double eps) {
        double[][] d2fXo={{8,-2},{-2,2}};
        double[][] invd2fXo={{1.0/6.0,1.0/6.0},{1.0/6.0,2.0/3.0}};
        while (absF(df(Xo)) > eps)
        {
            Xo[0]=Xo[0]-invd2fXo[0][0]*df(Xo)[0]-invd2fXo[0][1]*df(Xo)[1];
            Xo[1]=Xo[1]-invd2fXo[1][0]*df(Xo)[0]-invd2fXo[1][1]*df(Xo)[1];
        }
        return Xo;
    }

    public double[] prod(double[][] H, double[] f)
    {
        double[] res={H[0][0]*f[0]+H[0][1]*f[1],H[1][0]*f[0]+H[1][1]*f[1]};
        return res;
    }
    public double[] SR1Method(double[] Xo, double eps)
        {
        //---------------------------f(x)->min(alp) by Golden Section Method
        boolean flag = true;
        double a=-10,b=10,a2=a,b2=b,c,d,alp;
        double[] B={0,0},delta={0,0},gamma={0,0};
        double[][] H={{1,0},{0,1}};
        int iter=0,n=2;
        while(absF(df(Xo))>eps) {
            if(iter%n==0)
            {
                H[0][0]=1; H[0][1]=0; H[1][0]=0; H[1][1]=1;
            }
            else
            {
                double[][] h=H;
                H[0][0]=(delta[0]-H[0][0]*gamma[0]-H[0][1]*gamma[1])*(delta[0]-H[0][0]*gamma[0]-H[0][1]*gamma[1])/
                (gamma[0]*(delta[0]-H[0][0]*gamma[0]-H[0][1]*gamma[1])+(delta[1]-H[1][0]*gamma[0]-H[1][1]*gamma[1])*gamma[1]);
                H[0][1]=(delta[0]-H[0][0]*gamma[0]-H[0][1]*gamma[1])*(delta[1]-H[1][0]*gamma[0]-H[1][1]*gamma[1])/
                (gamma[0]*(delta[0]-H[0][0]*gamma[0]-H[0][1]*gamma[1])+(delta[1]-H[1][0]*gamma[0]-H[1][1]*gamma[1])*gamma[1]);
                H[1][0]=H[0][1];
                H[1][1]=(delta[1]-H[1][0]*gamma[0]-H[1][1]*gamma[1])*(delta[1]-H[1][0]*gamma[0]-H[1][1]*gamma[1])/
                (gamma[0]*(delta[0]-H[0][0]*gamma[0]-H[0][1]*gamma[1])+(delta[1]-H[1][0]*gamma[0]-H[1][1]*gamma[1])*gamma[1]);
            }
            flag=true;
            while (flag) {
                c = a + (b - a) * (3 - sqrt(5)) / 2;
                d = a + (b - a) * (sqrt(5) - 1) / 2;
                while (b - a >= eps*eps) {
                    B[0] = prod(H, df(Xo))[0];
                    B[1] = prod(H, df(Xo))[1];
                    if (Fxalp(Xo, B, c) <= Fxalp(Xo, B, d)) {
                        b = d;
                        d = c;
                        c = a + (b - a) * (3 - sqrt(5)) / 2;
                    } else {
                        a = c;
                        c = d;
                        d = a + (b - a) * (sqrt(5) - 1) / 2;
                    }
                }
                if ((a + b) / 2 >= a2 & (a + b) / 2 <= a2+eps) {
                    b = a + 2 * eps;
                    a = a - b2+a2;
                    a2=a;
                    b2=b;
                } else if ((a + b) / 2 >= b2-eps & (a + b) / 2 <= b2) {
                    a = b - 2 * eps;
                    b = b + b2-a2;
                    a2=a;
                    b2=b;
                } else {
                    flag = false;
                }
            }
            alp = -(a + b) / 2;
            double[] x = {0, 0};
            x[0] = Xo[0] - alp * prod(H, df(Xo))[0];
            x[1] = Xo[1] - alp * prod(H, df(Xo))[1];
            delta[0]=x[0]-Xo[0]; delta[1]=x[1]-Xo[1];
            gamma[0]=df(x)[0]-df(Xo)[0]; gamma[1]=df(x)[1]-df(Xo)[1];
            Xo = x;
            iter++;
        }
        return Xo;
    }
    public double[] ConjugateGradMethod(double[] Xo, double eps)
    {
        boolean flag = true;
        double a=-10,b=10,a2=a,b2=b,c,d,alp,beta=0;
        double[] D={-df(Xo)[0],-df(Xo)[1]},Xk={0,0};
        int iter=0,n=2;
        while(absF(df(Xo))>eps) {
            if(iter%n==0)
            {
                beta=0;
            }
            else
            {
                beta=(absF(df(Xo))*absF(df(Xo))) / (absF(df(Xk))*absF(df(Xk)));
            }
            D[0]=-df(Xo)[0]+beta*D[0];
            D[1]=-df(Xo)[1]+beta*D[1];
            flag=true;
            while (flag) {
                c = a + (b - a) * (3 - sqrt(5)) / 2;
                d = a + (b - a) * (sqrt(5) - 1) / 2;
                while (b - a >= eps*eps) {
                    if (Fxalp(Xo, D, c) <= Fxalp(Xo, D, d)) {
                        b = d;
                        d = c;
                        c = a + (b - a) * (3 - sqrt(5)) / 2;
                    } else {
                        a = c;
                        c = d;
                        d = a + (b - a) * (sqrt(5) - 1) / 2;
                    }
                }
                if ((a + b) / 2 >= a2 & (a + b) / 2 <= a2+eps) {
                    b = a + 2 * eps;
                    a = a - b2+a2;
                    a2=a;
                    b2=b;
                } else if ((a + b) / 2 >= b2-eps & (a + b) / 2 <= b2) {
                    a = b - 2 * eps;
                    b = b + b2-a2;
                    a2=a;
                    b2=b;
                } else {
                    flag = false;
                }
            }
            alp = (a + b) / 2;
            double[] x = {0, 0}; // x -  вспомогательная переменная
            x[0] = Xo[0] + alp * D[0];
            x[1] = Xo[1] + alp * D[1];
            Xk=Xo; // Xk - сохраняем значение точки с предыдущего шага для нахождения beta
            Xo = x;
            iter++;
        }
        return Xo;
    }
    public double f280(double[] x)
    {
        return 2*x[0]*x[0]+3*x[1]*x[1]-40*x[0]-48*x[1];
    }
    public double[] g(double[] x)
    {
        double[] res={-x[0], -x[1], x[0]-x[1]-6, x[1]+0.8*x[0]-12, x[1]-0.8*x[0]-4};
        return res;
    }
    public double H(double[] x)
    {
        double sum=0;
        for (int i=0;i<5;i++)
        {
            sum+= max(0,g(x)[i])*max(0,g(x)[i]);
        }
        return sum;
    }

    public double[] PenaltyMethod(double[] Xo, double eps)
    {
        //f(x)+r1*H(x) -> min
        // методом золотого сечения
        double a=-10, b=10, r=1;
        double[] Xc={0,0},Xd={0,0},x={Xo[0]+10,Xo[1]+10};
        do {
            do {
                x[0] = Xo[0];
                x[1] = Xo[1];
                for (int i = 0; i < 2; i++) {
                    a = -10;
                    b = 10;
                    double c = a + (b - a) * (3 - sqrt(5)) / 2, d = a + (b - a) * (sqrt(5) - 1) / 2;
                    while (b - a >= eps*eps) {
                        Xc[0] = Xo[0];
                        Xc[1] = Xo[1];
                        Xd[0] = Xo[0];
                        Xd[1] = Xo[1];
                        Xc[i] = Xo[i] + c;
                        Xd[i] = Xo[i] + d;
                        if (f280(Xc) + r * H(Xc) <= f280(Xd) + r * H(Xd)) {
                            b = d;
                            d = c;
                            c = a + (b - a) * (3 - sqrt(5)) / 2;
                        } else {
                            a = c;
                            c = d;
                            d = a + (b - a) * (sqrt(5) - 1) / 2;
                        }
                    }
                    Xo[i] += (a + b) / 2;
                }
            } while (abs(f280(x)+r*H(x) - f280(Xo)-r*H(Xo)) > eps*eps);
            r=r*10;
        } while(H(Xo)>eps);
        double y=f280(Xo);
        return Xo;
    }
    /** Called when the activity is first created. */
    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        // найдем View-элементы
        a = (TextView) findViewById(R.id.a);
        b = (TextView) findViewById(R.id.b);
        eps = (TextView) findViewById(R.id.eps);
        btnPassiveSearch = (Button) findViewById(R.id.btnPassiveSearch);
        btnDihotomySearch = (Button) findViewById(R.id.btnDihotomySearch);
        btnGoldenSection = (Button) findViewById(R.id.btnGoldenSection);
        btnFibonacciSearch = (Button) findViewById(R.id.btnFibonacciSearch);
        btnTangentMethod = (Button) findViewById(R.id.btnTangentMethod);
        btnNRMethod = (Button) findViewById(R.id.btnNRMethod);
        btnSecantMethod = (Button) findViewById(R.id.btnSecantMethod);
        enter_a = (EditText) findViewById(R.id.enter_a);
        enter_b = (EditText) findViewById(R.id.enter_b);
        enter_eps = (EditText) findViewById(R.id.enter_eps);
        xmin = (TextView) findViewById(R.id.xmin);
        ymin = (TextView) findViewById(R.id.ymin);
        X1 = (TextView) findViewById(R.id.X1);
        X2 = (TextView) findViewById(R.id.X2);
        L = (TextView) findViewById(R.id.L);
        alpha = (TextView) findViewById(R.id.alpha);
        enter_X1 = (EditText) findViewById(R.id.enter_X1);
        enter_X2 = (EditText) findViewById(R.id.enter_X2);
        enter_L = (EditText) findViewById(R.id.enter_L);
        enter_alpha = (EditText) findViewById(R.id.enter_alpha);
        btnSuccessiveDisplacementMethod = (Button) findViewById(R.id.btnSuccessiveDisplacementMethod);
        btnGradMethod= (Button) findViewById(R.id.btnGradMethod);
        btnGradConstStep= (Button) findViewById(R.id.btnGradConstStep);
        btnMNGS= (Button) findViewById(R.id.btnMNGS);
        btnDivergentSeriesMethod= (Button) findViewById(R.id.btnDivergentSeriesMethod);
        Xmin = (TextView) findViewById(R.id.Xmin);
        Ymin = (TextView) findViewById(R.id.Ymin);
        btnAcceleratedGradMethod = (Button) findViewById(R.id.btnAcceleratedGradMethod);
        btnRavineMethod = (Button) findViewById(R.id.btnRavineMethod);
        btnNR2Method = (Button) findViewById(R.id.btnNR2Method);
        btnConjugateGradMethod = (Button) findViewById(R.id.btnConjugateGradMethod);
        btnSR1Method = (Button) findViewById(R.id.btnSR1Method);
        Time = (TextView) findViewById(R.id.Time);
        btnPenaltyMethod = (Button) findViewById(R.id.btnPenaltyMethod);




        // создаем обработчик нажатия
        View.OnClickListener oclBtnPassiveSearch = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_a.getText().toString())
                        || TextUtils.isEmpty(enter_b.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double a1 = Float.parseFloat(enter_a.getText().toString());
                double b1 = Float.parseFloat(enter_b.getText().toString());
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double Xmin=PassiveSearch(a1, b1, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double Ymin = f(Xmin);
                xmin.setText(Double.toString(Xmin));
                ymin.setText(Double.toString(Ymin));
            }
        };
        View.OnClickListener oclBtnDihotomySearch = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_a.getText().toString())
                        || TextUtils.isEmpty(enter_b.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double a1 = Float.parseFloat(enter_a.getText().toString());
                double b1 = Float.parseFloat(enter_b.getText().toString());
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double Xmin= DihotomySearch(a1, b1, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double Ymin = f(Xmin);
                xmin.setText(Double.toString(Xmin));
                ymin.setText(Double.toString(Ymin));
            }
        };
        View.OnClickListener oclBtnGoldenSection = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_a.getText().toString())
                        || TextUtils.isEmpty(enter_b.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double a1 = Float.parseFloat(enter_a.getText().toString());
                double b1 = Float.parseFloat(enter_b.getText().toString());
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double Xmin= GoldenSection(a1, b1, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double Ymin = f(Xmin);
                xmin.setText(Double.toString(Xmin));
                ymin.setText(Double.toString(Ymin));
            }
        };
        View.OnClickListener oclBtnFibonacciSearch = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_a.getText().toString())
                        || TextUtils.isEmpty(enter_b.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double a1 = Float.parseFloat(enter_a.getText().toString());
                double b1 = Float.parseFloat(enter_b.getText().toString());
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double Xmin= FibonacciSearch(a1, b1, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double Ymin = f(Xmin);
                xmin.setText(Double.toString(Xmin));
                ymin.setText(Double.toString(Ymin));
            }
        };
        View.OnClickListener oclBtnTangentMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_a.getText().toString())
                        || TextUtils.isEmpty(enter_b.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double a1 = Float.parseFloat(enter_a.getText().toString());
                double b1 = Float.parseFloat(enter_b.getText().toString());
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double Xmin= TangentMethod(a1, b1, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double Ymin = f(Xmin);

                xmin.setText(Double.toString(Xmin));
                ymin.setText(Double.toString(Ymin));
            }
        };
        View.OnClickListener oclBtnNewtonRaphsonMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_a.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double a1 = Float.parseFloat(enter_a.getText().toString());
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double Xmin= NewtonRaphsonMethod(a1, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double Ymin = f(Xmin);
                xmin.setText(Double.toString(Xmin));
                ymin.setText(Double.toString(Ymin));
            }
        };
        View.OnClickListener oclBtnSecantMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_a.getText().toString())
                        || TextUtils.isEmpty(enter_b.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double a1 = Float.parseFloat(enter_a.getText().toString());
                double b1 = Float.parseFloat(enter_b.getText().toString());
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double Xmin= SecantMethod(a1, b1, eps), Ymin = f(Xmin);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                xmin.setText(Double.toString(Xmin));
                ymin.setText(Double.toString(Ymin));
            }
        };
        View.OnClickListener oclBtnSuccessiveDisplacementMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double[] xmin=SuccessiveDisplacementMethod(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnGradMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_L.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                double L = Float.parseFloat(enter_L.getText().toString());
                System.out.printf("BEFORE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%n");
                long startTime = System.currentTimeMillis();
                double[] xmin=GradMethod(x, eps, L);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnGradConstStep = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_alpha.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                double alpha = Float.parseFloat(enter_alpha.getText().toString());
                System.out.printf("BEFORE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%n");
                long startTime = System.currentTimeMillis();
                double[] xmin=GradConstStep(x, eps, alpha);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnMNGS = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                System.out.printf("BEFORE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%n");
                long startTime = System.currentTimeMillis();
                double[] xmin=MNGS(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnDivergentSeriesMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                System.out.printf("BEFORE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%n");
                long startTime = System.currentTimeMillis();
                double[] xmin=DivergentSeriesMethod(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnAcceleratedGradMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                System.out.printf("BEFORE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%n");
                long startTime = System.currentTimeMillis();
                double[] xmin=AcceleratedGradMethod(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnRavineMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                System.out.printf("BEFORE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%n");
                long startTime = System.currentTimeMillis();
                double[] xmin=RavineMethod(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnNR2Method = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double[] xmin=NR2Method(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnSR1Method = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double[] xmin=SR1Method(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));
            }
        };
        View.OnClickListener oclBtnConjugateGradMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double[] xmin=ConjugateGradMethod(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));

            }
        };
        View.OnClickListener oclBtnPenaltyMethod = new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                // Проверяем поля на пустоту
                if (TextUtils.isEmpty(enter_X1.getText().toString())
                        || TextUtils.isEmpty(enter_X2.getText().toString())
                        || TextUtils.isEmpty(enter_eps.getText().toString())) {
                    return;
                }
                // читаем EditText и заполняем переменные числами
                double x1 = Float.parseFloat(enter_X1.getText().toString());
                double x2 = Float.parseFloat(enter_X2.getText().toString());
                double[] x={x1,x2};
                double eps = Float.parseFloat(enter_eps.getText().toString());
                long startTime = System.currentTimeMillis();
                double[] xmin=PenaltyMethod(x, eps);
                long timeSpent = System.currentTimeMillis() - startTime;
                Time.setText(Double.toString(timeSpent));
                System.out.println("программа выполнялась " + timeSpent + " милисекунд");
                double ymin = f280(xmin);
                Xmin.setText(Double.toString(xmin[0])+"; "+Double.toString(xmin[1]));
                Ymin.setText(Double.toString(ymin));

            }
        };
        btnPassiveSearch.setOnClickListener(oclBtnPassiveSearch);
        btnDihotomySearch.setOnClickListener(oclBtnDihotomySearch);
        btnGoldenSection.setOnClickListener(oclBtnGoldenSection);
        btnFibonacciSearch.setOnClickListener(oclBtnFibonacciSearch);
        btnTangentMethod.setOnClickListener(oclBtnTangentMethod);
        btnNRMethod.setOnClickListener(oclBtnNewtonRaphsonMethod);
        btnSecantMethod.setOnClickListener(oclBtnSecantMethod);
        btnSuccessiveDisplacementMethod.setOnClickListener(oclBtnSuccessiveDisplacementMethod);
        btnGradMethod.setOnClickListener(oclBtnGradMethod);
        btnGradConstStep.setOnClickListener(oclBtnGradConstStep);
        btnMNGS.setOnClickListener(oclBtnMNGS);
        btnDivergentSeriesMethod.setOnClickListener(oclBtnDivergentSeriesMethod);
        btnAcceleratedGradMethod.setOnClickListener(oclBtnAcceleratedGradMethod);
        btnRavineMethod.setOnClickListener(oclBtnRavineMethod);
        btnNR2Method.setOnClickListener(oclBtnNR2Method);
        btnSR1Method.setOnClickListener(oclBtnSR1Method);
        btnConjugateGradMethod.setOnClickListener(oclBtnConjugateGradMethod);
        btnPenaltyMethod.setOnClickListener(oclBtnPenaltyMethod);
    }
}
