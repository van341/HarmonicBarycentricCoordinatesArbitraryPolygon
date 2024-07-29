
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Numerics;
using VectorDouble = MathNet.Numerics.LinearAlgebra.Vector<double>;
using MathNet.Numerics;

namespace barycentric_calc
{

    internal sealed partial class LegenderePolynomial
    {
        private int N;
        private VectorDouble coefficient1;
        private VectorDouble coefficient2;
        private Dictionary<int, int> factorialCache = new Dictionary<int, int>();

        private int Factorial(int n)
        {
            if (factorialCache.ContainsKey(n))
            {
                return factorialCache[n];
            }

            int f = 1;
            for (int i = 2; i <= n; ++i)
            {
                f *= i;
            }

            factorialCache[n] = f;
            return f;
        }

        public LegenderePolynomial(int N)
        {
            this.N = N;
            coefficient1 = new DenseVector(N);
            coefficient2 = new DenseVector(N);
            for (int n = 0; n < N; ++n)
            {
                coefficient1[n] = (2 * n + 1.0) / (n + 1.0);
                coefficient2[n] = n / (n + 1.0);
            }
        }

        public double FirstKindD(int n, double x)
        {
            double p0 = 1.0;
            if (n == 0)
            {
                return p0;
            }
            int i = 1;
            double p1 = x;
            while (i < n)
            {
                (p0, p1) = (p1, p0);
                p1 = coefficient1[i] * x * p0 - coefficient2[i] * p1;
                ++i;
            }
            return p1;
        }

        public Complex FirstKind(int n, Complex x)
        {
            Complex p0 = Complex.One;
            if (n == 0)
            {
                return p0;
            }
            int i = 1;
            Complex p1 = x;
            while (i < n)
            {
                (p0, p1) = (p1, p0);
                p1 = coefficient1[i] * x * p0 - coefficient2[i] * p1;
                ++i;
            }
            return p1;
        }

        public Complex SecondKind(int n, Complex x)
        {
            Complex p0 = 0.5 * Complex.Log((x + 1.0) / (x - 1.0));
            if (n == 0)
            {
                return p0;
            }
            int i = 1;
            Complex p1 = x * p0 - 1.0;
            while (i < n)
            {
                (p0, p1) = (p1, p0);
                p1 = coefficient1[i] * x * p0 - coefficient2[i] * p1;
                ++i;
            }
            return p1;
        }

        public double SecondKindD(int n, double x)
        {
            double p0 = 0.5 * Math.Log((x + 1.0) / (1.0 - x));
            if (n == 0)
            {
                return p0;
            }
            int i = 1;
            double p1 = x * p0 - 1.0;
            while (i < n)
            {
                (p0, p1) = (p1, p0);
                p1 = coefficient1[i] * x * p0 - coefficient2[i] * p1;
                ++i;
            }
            return p1;
        }

        public List<Complex> SecondKindVec(int n, Complex x)
        {
            var v = new List<Complex>();
            Complex p0 = 0.5 * Complex.Log((x + 1.0) / (x - 1.0));
            v.Add(p0);
            if (n == 0)
            {
                return v;
            }
            int i = 1;
            Complex p1 = x * p0 - 1.0;
            v.Add(p1);
            while (i < n)
            {
                (p0, p1) = (p1, p0);
                p1 = coefficient1[i] * x * p0 - coefficient2[i] * p1;
                ++i;
                v.Add(p1);
            }
            return v;
        }

        public List<Complex> SecondKindVec2(int n, Complex x)
        {
            var v = new List<Complex>();
            Complex p0 = -Complex.Log((x + 1.0) / (x - 1.0));
            v.Add(p0);
            if (n == 0)
            {
                return v;
            }
            int i = 1;
            Complex p1 = x * p0 + 2.0;
            v.Add(p1);
            while (i < n)
            {
                (p0, p1) = (p1, p0);
                p1 = coefficient1[i] * x * p0 - coefficient2[i] * p1;
                ++i;
                v.Add(p1);
            }
            return v;
        }

        public double DerivativeFirstKind(int n, int m, double x)
        {
            if (n < m)
            {
                return 0;
            }
            if (n == m)
            {
                return Factorial(2 * n) / (Factorial(n) * Math.Pow(2, n));
            }
            if (n - 1 == m)
            {
                return x * Factorial(2 * n) / (Factorial(n) * Math.Pow(2, n));
            }

            double p0 = 1.0 / (Factorial(m) * Math.Pow(2, m)) * Factorial(2 * m);
            double p1 = 1.0 / (Factorial(m + 1) * Math.Pow(2, m + 1)) * x * Factorial(2 * (m + 1));
            int i = m + 1;
            while (i < n)
            {
                (p0, p1) = (p1, p0);
                p1 = ((2.0 * i + 1.0) * x * p0 - (i + m) * p1) / (i + 1.0 - m);
                ++i;
            }
            return p1;
        }
    }

    internal sealed partial class computation_zeros_legendre_polynomials
    {
        private int N;
        private int K;
        private LegenderePolynomial legendere;

        public VectorDouble x0;
        public VectorDouble A0;

        public computation_zeros_legendre_polynomials(int N, int K)
        {
            this.N = N;
            this.K = K;
            legendere = new LegenderePolynomial(N);
            x0 = new DenseVector(N);
            A0 = new DenseVector(N);
            for (int n = 0; n < N; ++n)
            {
                x0[n] = Math.Cos((Math.PI * (4.0 * (n + 1) - 1)) / (4.0 * N + 2));
            }
            x0 = iterative_refinement_Newton_method(x0);
            double tr;
            for (int n = 0; n < N; ++n)
            {
                tr = legendere.DerivativeFirstKind(N, 1, x0[n]);
                A0[n] = 2.0 / ((1 - x0[n] * x0[n]) * tr * tr);
            }
        }

        private VectorDouble iterative_refinement_Newton_method(VectorDouble x0)
        {
            var x = new DenseVector(N);
            x0.CopyTo(x);
            for (int k = 0; k < K; ++k)
            {
                for (int n = 0; n < N; ++n)
                {
                    double firstKindD = legendere.FirstKindD(N, x[n]);
                    double derivativeFirstKind = legendere.DerivativeFirstKind(N, 1, x[n]);
                    x[n] = x[n] - firstKindD / derivativeFirstKind;
                }
            }
            return x;
        }

    }

    internal sealed partial class barycentric_calculations_interier {
        private List<Complex> P;
        private List<Complex> L;
        private List<Matrix<double>> Coeff;
        private int NK;
        private int MN;
	    private int MM;
        private LegenderePolynomial legendere;
        private computation_zeros_legendre_polynomials zeros_legendre;
        public barycentric_calculations_interier(List<Complex> P, int NK, int MN, int MM)
        {
            this.P = new List<Complex>(P);
            this.NK = NK;
            this.MN = MN;
            this.MM = MM;
            L = calculations_L(P);
            legendere = new LegenderePolynomial(NK);
            var CC = weight_coefficients(NK, MM, MN);
            int N = P.Count;
            int NL = NK * N;
            Coeff = new List<Matrix<double>>();

            // Предварительные вычисления для ii1 и hn
            int[] ii1Array = new int[NL];
            int[] hnArray = new int[NL];
            for (int n = 0; n < NL; ++n)
            {
                ii1Array[n] = (int)Math.Ceiling((n + 1.0) / NK) - 1;
                hnArray[n] = sp(n, NK - 1);
            }

            for (int j = 0; j < N; ++j)
            {
                Coeff.Add(new DenseMatrix(N, NK));
                for (int n = 0; n < NL; ++n)
                {
                    Coeff[j][ii1Array[n], hnArray[n]] = CC[j][n];
                }
            }
        }

        private int sp(int i, int n)
        {
            return (i >= 0) ? (i % (n + 1)) : (n + i + 1);
        }

        private List<Complex> calculations_L(List<Complex> P)
        {
            int N = P.Count;
            var LL = new List<Complex>();
            for (int i = 0; i < N; ++i)
            {
                int j = sp(i + 1, N - 1);
                LL.Add(P[j] - P[i]);
            }
            return LL;
        }

        private double function(int n, int m, Complex W, Complex D, double t)
	    {
		    Complex ff = W * (t + 1.0) + D;
            double result = legendere.SecondKind(m, ff).Imaginary;
            result *= legendere.FirstKindD(n, t);
		    return result;
	    }

        private double int_eigenvaluecorel(int i, int j, int n, int m, Matrix<Complex> W, Matrix<Complex> D, Matrix<double> X, VectorDouble A) {
		double result = 0;
		if (i != j) {
			for (int t1 = 0; t1<X.RowCount; ++t1) {
				double sum = 0;
				for (int m1 = 0; m1<X.ColumnCount; ++m1) {
					sum += function(n, m, W[i, j], D[i, j], X[t1, m1]) * A[m1];
                }
                result += sum;
			}
            result /= -Math.PI;
            result /= X.RowCount;
		    }
		    return result;
	    }

        private double boundary_conditions(int k, int j, int n)
        {
            double result = 0;
            int N = P.Count;
            if (k == j)
            {
                if (n == 0) { result = 1.0; }
                else
                {
                    if (n == 1)
                    {
                        result = Math.Sqrt(3.0);
                        result = -1.0 / result;
                    }
                }
            }
            else
            {
                int j1 = sp(j - 2, N - 1) + 1;
                if (j1 == k)
                {
                    if (n == 0) { result = 1.0; }
                    else
                    {
                        if (n == 1)
                        {
                            result = Math.Sqrt(3.0);
                            result = 1.0 / result;
                        }
                    }
                }
            }
            return result;
        }

        private List<VectorDouble> weight_coefficients(int NK, int N0_L, int MN) 
        {
		    int N = P.Count;
            int NL = N * NK;
            int ii1, jj1, hn, hm;
            var A = new DenseMatrix(NL, NL);
            int K0_L = 4000;
            zeros_legendre = new computation_zeros_legendre_polynomials(N0_L, K0_L);
		    var X_ = new DenseMatrix(MN, N0_L);
            double deltam = 1.0 / MN;
		    for (int n = 0; n<MN; ++n) {
			    for (int k = 0; k<N0_L; ++k) {
				    X_[n, k] = 2.0* deltam* (0.5 * (zeros_legendre.x0[k] + 1.0) + n) - 1.0;
			    }
            }
            
            var W = new MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix(L.Count, L.Count);
            var D = new MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix(L.Count, L.Count);
            for (int n = 0; n < L.Count; ++n)
            {
                for (int k = 0; k < L.Count; ++k)
                {
                    if (n != k)
                    {
                        W[n, k] = L[n] / L[k];
                        D[n, k] = 2.0 * (P[n] - P[k]) / L[k] - 1.0;
                    }
                    else
                    {
                        W[n, k] = 1;
                        D[n, k] = 1;
                    }
                }
            }
            var wk = new DenseVector(NK);
            for (int k = 0; k < NK; ++k)
            {
                wk[k] = Math.Sqrt(2.0 * k + 1.0);
            }

            //расчет в один поток
            for (int n = 0; n < NL; ++n)
            {
                ii1 = (int)Math.Ceiling((n + 1.0) / NK) - 1;
                hn = sp(n, NK - 1);
                for (int m = 0; m < NL; ++m)
                {
                    jj1 = (int)Math.Ceiling((m + 1.0) / NK) - 1;
                    hm = sp(m, NK - 1);
                    if (n != m)
                    {
                        A[n, m] = int_eigenvaluecorel(ii1, jj1, hn, hm, W, D, X_, zeros_legendre.A0);
                        A[n, m] *= wk[hn] * wk[hm];
                    }
                    else
                    {
                        A[n, m] = 1;
                    }
                }
            }

            var C = new List<VectorDouble>();
            var F = new List<VectorDouble>();
            Matrix<double> AC = A.Inverse();
            for (int k = 0; k < N; ++k)
            {
                F.Add(new DenseVector(NL));
                C.Add(new DenseVector(NL));
                for (int n = 0; n < NL; ++n)
                {
                    ii1 = (int)Math.Ceiling((n + 1.0) / NK);
                    hn = sp(n, NK - 1);
                    F[k][n] = boundary_conditions(ii1, k + 1, hn);
                }
                C[k] = AC * F[k];
            }
            return C;
	    }
     

        private List<double> lambda_vec(Complex z, int n)
        {
            var zz = legendere.SecondKindVec2(n, z);
            var R = new List<double>();
            for (int i = 0; i < n; ++i)
            {
                R.Add(Math.Sqrt(2 * i + 1.0) * zz[i].Imaginary);
            }
            return R;
        }

        private void dlambda_xy(Complex z, int n, int i, List<double> R, List<double> dR_x, List<double> dR_y)
        {
            var leg = legendere.SecondKindVec(n, z);
            var zz = new List<Complex>();
            Complex rr = 1.0 / (1.0 - z * z);
            for (int j = 0; j < n; ++j)
            {
                if (j == 0)
                {
                    zz.Add(rr / L[i]);
                }
                else
                {
                    zz.Add(j*(rr * (leg[j - 1] - z * leg[j]) / L[i]));
                }
            }
            double tt;
            for (int j = 0; j < n; ++j)
            {
                tt = -2.0* Math.Sqrt(2 * j + 1.0);
                R.Add(tt * leg[j].Imaginary);
                dR_x.Add(2 * tt * zz[j].Imaginary);
                dR_y.Add(2 * tt * zz[j].Real);
            }
        }

        public List<double> barycentric(Complex z)
        {
            int N = P.Count;
            var LL_ = new DenseMatrix(N, NK);
            var zz = new List<Complex>();
            for (int k = 0; k < N; ++k)
            {
                zz.Add(2.0* (z - P[k]) / L[k] - 1.0);
            }
            for (int k = 0; k < N; ++k)
            {
                var ss = lambda_vec(zz[k], NK);
                for (int n = 0; n < NK; ++n)
                {
                    LL_[k, n] = ss[n];
                }
            }
            var zeta = new List<double>();
            var ss_ = new DenseVector(N);
            for (int j = 0; j < N; ++j)
            {
                for (int k = 0; k < N; ++k)
                {
                    ss_[k] = Coeff[j].Row(k).DotProduct(LL_.Row(k));
                }
                zeta.Add(ss_.Sum() / (2.0 * Math.PI));
            }
            return zeta;
        }

        public void dxy_barycentric(Complex z, List<double> b, List<double> db_x, List<double> db_y)
        {
            int N = P.Count;
            var LL_ = new DenseMatrix(N, NK);
            var LL_x = new DenseMatrix(N, NK);
            var LL_y = new DenseMatrix(N, NK);
            var zz = new List<Complex>();
            for (int k = 0; k < N; ++k)
            {
                zz.Add(2.0* (z - P[k]) / L[k] - 1.0);
            }
            var R = new List<double>();
            var dR_x = new List<double>();
            var dR_y = new List<double>();
            for (int k = 0; k < N; ++k)
            {
                dlambda_xy(zz[k], NK, k, R, dR_x, dR_y);
                for (int n = 0; n < NK; ++n)
                {
                    LL_[k, n] = R[n];
                    LL_x[k, n] = dR_x[n];
                    LL_y[k, n] = dR_y[n];
                }
            }
            var ss = new DenseVector(N);
            var ssx = new DenseVector(N);
            var ssy = new DenseVector(N);   
            for (int j = 0; j < N; ++j)
            {
                for (int k = 0; k < N; ++k)
                {
                    ss[k] = Coeff[j].Row(k).DotProduct(LL_.Row(k));
                    ssx[k] = Coeff[j].Row(k).DotProduct(LL_x.Row(k));
                    ssy[k] = Coeff[j].Row(k).DotProduct(LL_y.Row(k));
                }
                b.Add(ss.Sum() / (2.0* Math.PI));
                db_x.Add(ssx.Sum() / (2.0 * Math.PI));
                db_y.Add(ssy.Sum() / (2.0 * Math.PI));
            }
        }

        public bool check_point_in_polygon(Complex p)
	    {
		    bool result = false;
		    for (int i1 = 0; i1<P.Count; ++i1)
		    {
			    int i2 = (i1 + 1) % P.Count;
			    if (((p.Imaginary > P[i1].Imaginary) ^ (p.Imaginary > P[i2].Imaginary)) &&
				    (p.Real > P[i1].Real + (L[i1].Real)*(p.Imaginary - P[i1].Imaginary) / (L[i1].Imaginary)))
				    result = !result;
		    }
		    return result;
	    }

    }

    internal sealed partial class barycentric_calculations_exterier
    {
        private List<Complex> P;
        private List<Complex> L;
        private List<double> Length_L;
        private List<Matrix<double>> Coeff;
        private int NK;
        private int MN;
        private int MM;
        private LegenderePolynomial legendere;
        private computation_zeros_legendre_polynomials zeros_legendre;
        public barycentric_calculations_exterier(List<Complex> P, int NK, int MN, int MM)
        {
            this.P = new List<Complex>(P);
            this.NK = NK;
            this.MN = MN;
            this.MM = MM;
            L = calculations_L(P);
            Length_L = new List<double>();
            for (int i = 0; i < L.Count; ++i)
            {
                Length_L.Add(L[i].Norm());
            }
            legendere = new LegenderePolynomial(NK);
            var CC = weight_coefficients(NK, MM, MN);
            int N = P.Count;
            int NL = NK * N;
            Coeff = new List<Matrix<double>>();
            int ii1, hn;
            for (int j = 0; j < N; ++j)
            {
                Coeff.Add(new DenseMatrix(N, NK));
                for (int n = 0; n < NL; ++n)
                {
                    ii1 = (int)Math.Ceiling((n + 1.0) / NK) - 1;
                    hn = sp(n, NK - 1);
                    Coeff[j][ii1, hn] = CC[j][n];
                }
            }
        }

        private int sp(int i, int n)
        {
            int m;
            if (i >= 0)
            {
                m = (i) % (n + 1);
            }
            else
            {
                m = n + i + 1;
            }
            return m;
        }

        private List<Complex> calculations_L(List<Complex> P)
        {
            int N = P.Count;
            var LL = new List<Complex>();
            for (int i = 0; i < N; ++i)
            {
                int j = sp(i + 1, N - 1);
                LL.Add(P[j] - P[i]);
            }
            return LL;
        }

        private double function(int n, int m, Complex W, Complex D, double t)
        {
            Complex ff = W * (t + 1.0) + D;
            double result = legendere.SecondKind(m, ff).Imaginary;
            result *= legendere.FirstKindD(n, t);
            return result;
        }

        private double int_eigenvaluecorel(int i, int j, int n, int m, Matrix<Complex> W, Matrix<Complex> D, Matrix<double> X, VectorDouble A)
        {
            double result = 0;
            if (i != j)
            {
                for (int t1 = 0; t1 < X.RowCount; ++t1)
                {
                    double sum = 0;
                    for (int m1 = 0; m1 < X.ColumnCount; ++m1)
                    {
                        sum += function(n, m, W[i, j], D[i, j], X[t1, m1]) * A[m1];
                    }
                    result += sum;
                }
                result /= -Math.PI;
                result /= X.RowCount;
            }
            return result;
        }

        private double boundary_conditions(int k, int j, int n)
        {
            double result = 0;
            int N = P.Count;
            if (k == j)
            {
                if (n == 0) { result = 1.0; }
                else
                {
                    if (n == 1)
                    {
                        result = 3.0;
                        result = Math.Sqrt(result);
                        result = -1.0 / result;
                    }
                }
            }
            else
            {
                int j1 = sp(j - 2, N - 1) + 1;
                if (j1 == k)
                {
                    if (n == 0) { result = 1.0; }
                    else
                    {
                        if (n == 1)
                        {
                            result = 3.0;
                            result = Math.Sqrt(result);
                            result = 1.0 / result;
                        }
                    }
                }
            }
            return result;
        }

        private List<VectorDouble> weight_coefficients(int NK, int N0_L, int MN)
        {
            int N = P.Count;
            int NL = N * NK;
            int ii1, jj1, hn, hm;
            var A = new DenseMatrix(NL, NL);
            int K0_L = 4000;
            zeros_legendre = new computation_zeros_legendre_polynomials(N0_L, K0_L);
            var X_ = new DenseMatrix(MN, N0_L);
            double deltam = 1.0 / MN;
            for (int n = 0; n < MN; ++n)
            {
                for (int k = 0; k < N0_L; ++k)
                {
                    X_[n, k] = 2.0 * deltam * (0.5 * (zeros_legendre.x0[k] + 1.0) + n) - 1.0;
                }
            }

            var W = new MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix(L.Count, L.Count);
            var D = new MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix(L.Count, L.Count);
            for (int n = 0; n < L.Count; ++n)
            {
                for (int k = 0; k < L.Count; ++k)
                {
                    if (n != k)
                    {
                        W[n, k] = L[n] / L[k];
                        D[n, k] = 2.0 * (P[n] - P[k]) / L[k] - 1.0;
                    }
                    else
                    {
                        W[n, k] = 1;
                        D[n, k] = 1;
                    }
                }
            }
            var wk = new DenseVector(NK);
            for (int k = 0; k < NK; ++k)
            {
                wk[k] = Math.Sqrt(2.0 * k + 1.0);
            }

            //расчет в один поток
            for (int n = 0; n < NL; ++n)
            {
                ii1 = (int)Math.Ceiling((n + 1.0) / NK) - 1;
                hn = sp(n, NK - 1);
                for (int m = 0; m < NL; ++m)
                {
                    jj1 = (int)Math.Ceiling((m + 1.0) / NK) - 1;
                    hm = sp(m, NK - 1);
                    if (n != m)
                    {
                        A[n, m] = -int_eigenvaluecorel(ii1, jj1, hn, hm, W, D, X_, zeros_legendre.A0);
                        A[n, m] *= wk[hn] * wk[hm];
                    }
                    else
                    {
                        A[n, m] = 1;
                    }
                    if ((hn == 0) && (hm == 0)) { A[n, m] -= 2.0 * Length_L[jj1]; }
                }
            }

            var C = new List<VectorDouble>();
            var F = new List<VectorDouble>();
            Matrix<double> AC = A.Inverse();
            for (int k = 0; k < N; ++k)
            {
                F.Add(new DenseVector(NL));
                C.Add(new DenseVector(NL));
                for (int n = 0; n < NL; ++n)
                {
                    ii1 = (int)Math.Ceiling((n + 1.0) / NK);
                    hn = sp(n, NK - 1);
                    F[k][n] = boundary_conditions(ii1, k + 1, hn);
                }
                C[k] = -AC * F[k];
            }
            return C;
        }

        private List<double> lambda_vec(Complex z, int n)
        {
            var zz = legendere.SecondKindVec2(n, z);
            var R = new List<double>();
            for (int i = 0; i < n; ++i)
            {
                R.Add(Math.Sqrt(2 * i + 1.0) * zz[i].Imaginary);
            }
            return R;
        }

        private void dlambda_xy(Complex z, int n, int i, List<double> R, List<double> dR_x, List<double> dR_y)
        {
            var leg = legendere.SecondKindVec(n, z);
            var zz = new List<Complex>();
            Complex rr = 1.0 / (1.0 - z * z);
            for (int j = 0; j < n; ++j)
            {
                if (j == 0)
                {
                    zz.Add(rr / L[i]);
                }
                else
                {
                    zz.Add(j * (rr * (leg[j - 1] - z * leg[j]) / L[i]));
                }
            }
            double tt;
            for (int j = 0; j < n; ++j)
            {
                tt = -2.0 * Math.Sqrt(2 * j + 1.0);
                R.Add(tt * leg[j].Imaginary);
                dR_x.Add(2 * tt * zz[j].Imaginary);
                dR_y.Add(2 * tt * zz[j].Real);
            }
        }

        public List<double> barycentric(Complex z)
        {
            int N = P.Count;
            var LL_ = new DenseMatrix(N, NK);
            var zz = new List<Complex>();
            for (int k = 0; k < N; ++k)
            {
                zz.Add(2.0 * (z - P[k]) / L[k] - 1.0);
            }
            for (int k = 0; k < N; ++k)
            {
                var ss = lambda_vec(zz[k], NK);
                ss[0] += 2.0 * Math.PI * Length_L[k];
                for (int n = 0; n < NK; ++n)
                {
                    LL_[k, n] = ss[n];
                }
            }
            var zeta = new List<double>();
            var ss_ = new DenseVector(N);
            for (int j = 0; j < N; ++j)
            {
                for (int k = 0; k < N; ++k)
                {
                    ss_[k] = Coeff[j].Row(k).DotProduct(LL_.Row(k));
                }
                zeta.Add(ss_.Sum() / (2.0 * Math.PI));
            }
            return zeta;
        }

        public void dxy_barycentric(Complex z, List<double> b, List<double> db_x, List<double> db_y)
        {
            int N = P.Count;
            var LL_ = new DenseMatrix(N, NK);
            var LL_x = new DenseMatrix(N, NK);
            var LL_y = new DenseMatrix(N, NK);
            var zz = new List<Complex>();
            for (int k = 0; k < N; ++k)
            {
                zz.Add(2.0 * (z - P[k]) / L[k] - 1.0);
            }
            var R = new List<double>();
            var dR_x = new List<double>();
            var dR_y = new List<double>();
            for (int k = 0; k < N; ++k)
            {
                dlambda_xy(zz[k], NK, k, R, dR_x, dR_y);
                for (int n = 0; n < NK; ++n)
                {
                    LL_[k, n] = R[n];
                    LL_x[k, n] = dR_x[n];
                    LL_y[k, n] = dR_y[n];
                }
            }
            var ss = new DenseVector(N);
            var ssx = new DenseVector(N);
            var ssy = new DenseVector(N);
            for (int j = 0; j < N; ++j)
            {
                for (int k = 0; k < N; ++k)
                {
                    ss[k] = Coeff[j].Row(k).DotProduct(LL_.Row(k));
                    ssx[k] = Coeff[j].Row(k).DotProduct(LL_x.Row(k));
                    ssy[k] = Coeff[j].Row(k).DotProduct(LL_y.Row(k));
                }
                b.Add(ss.Sum() / (2.0 * Math.PI));
                db_x.Add(ssx.Sum() / (2.0 * Math.PI));
                db_y.Add(ssy.Sum() / (2.0 * Math.PI));
            }
        }

    }

    internal sealed partial class barycentric_generalized_calculations
    {
        private List<Complex> P;
        private List<Complex> L;
        private List<double> Length_L;
        private List<Matrix<double>> Coeff_i;
        private List<Matrix<double>> Coeff_e;
        private int NK;
        private int MN;
        private int MM;
        private LegenderePolynomial legendere;
        private computation_zeros_legendre_polynomials zeros_legendre;

        public barycentric_generalized_calculations(List<Complex> P, int NK, int MN, int MM)
        {
            this.P = new List<Complex>(P);
            this.NK = NK;
            this.MN = MN;
            this.MM = MM;
            L = calculations_L(P);
            Length_L = new List<double>(L.Count);
            for (int i = 0; i < L.Count; ++i)
            {
                Length_L.Add(L[i].Norm());
            }
            legendere = new LegenderePolynomial(NK);
            var CC = weight_coefficients_interier_and_exterier(NK, MM, MN);
            int N = P.Count;
            int NL = NK * N;
            Coeff_i = new List<Matrix<double>>(N);
            Coeff_e = new List<Matrix<double>>(N);

            for (int j = 0; j < N; ++j)
            {
                var coeff_i_j = new DenseMatrix(N, NK);
                var coeff_e_j = new DenseMatrix(N, NK);
                for (int n = 0; n < NL; ++n)
                {
                    int ii1 = (n / NK);
                    int hn = sp(n, NK - 1);
                    coeff_i_j[ii1, hn] = CC.Ci_[j][n];
                    coeff_e_j[ii1, hn] = CC.Ce_[j][n];
                }
                Coeff_i.Add(coeff_i_j);
                Coeff_e.Add(coeff_e_j);
            }
        }

        private struct Coef
        {
            public List<VectorDouble> Ci_;
            public List<VectorDouble> Ce_;

        };

        private int sp(int i, int n)
        {
            int m;
            if (i >= 0)
            {
                m = (i) % (n + 1);
            }
            else
            {
                m = n + i + 1;
            }
            return m;
        }


        private List<Complex> calculations_L(List<Complex> P)
        {
            int N = P.Count;
            var LL = new List<Complex>();
            for (int i = 0; i < N; ++i)
            {
                int j = sp(i + 1, N - 1);
                LL.Add(P[j] - P[i]);
            }
            return LL;
        }

        private double function(int n, int m, Complex W, Complex D, double t)
        {
            Complex ff = W * (t + 1.0) + D;
            double result = legendere.SecondKind(m, ff).Imaginary;
            result *= legendere.FirstKindD(n, t);
            return result;
        }

        private double int_eigenvaluecorel(int i, int j, int n, int m, Matrix<Complex> W, Matrix<Complex> D, Matrix<double> X, VectorDouble A)
        {
            double result = 0;
            if (i != j)
            {
                for (int t1 = 0; t1 < X.RowCount; ++t1)
                {
                    double sum = 0;
                    for (int m1 = 0; m1 < X.ColumnCount; ++m1)
                    {
                        sum += function(n, m, W[i, j], D[i, j], X[t1, m1]) * A[m1];
                    }
                    result += sum;
                }
                result /= -Math.PI;
                result /= X.RowCount;
            }
            return result;
        }
        private double boundary_conditions(int k, int j, int n)
        {
            if (n != 0 && n != 1)
            {
                return 0; // Обработка случаев, когда n не равно 0 или 1
            }

            if (k == j || sp(j - 2, P.Count - 1) + 1 == k)
            {
                if (n == 0)
                {
                    return 1.0;
                }
                else // n == 1
                {
                    return (k == j) ? -1.0 / Math.Sqrt(3.0) : 1.0 / Math.Sqrt(3.0);
                }
            }

            return 0;
        }

        private Coef weight_coefficients_interier_and_exterier(int NK, int N0_L, int MN)
        {
            int N = P.Count;
            int NL = N * NK;
            int ii1, jj1, hn, hm;
            var Ai = new DenseMatrix(NL, NL);
            var Ae = new DenseMatrix(NL, NL);
            int K0_L = 600;
            zeros_legendre = new computation_zeros_legendre_polynomials(N0_L, K0_L);
            var X_ = new DenseMatrix(MN, N0_L);
            double deltam = 1.0 / MN;
            for (int n = 0; n < MN; ++n)
            {
                for (int k = 0; k < N0_L; ++k)
                {
                    X_[n, k] = 2.0 * deltam * (0.5 * (zeros_legendre.x0[k] + 1.0) + n) - 1.0;
                }
            }

            var W = new MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix(L.Count, L.Count);
            var D = new MathNet.Numerics.LinearAlgebra.Complex.DenseMatrix(L.Count, L.Count);
            for (int n = 0; n < L.Count; ++n)
            {
                for (int k = 0; k < L.Count; ++k)
                {
                    if (n != k)
                    {
                        W[n, k] = L[n] / L[k];
                        D[n, k] = 2.0 * (P[n] - P[k]) / L[k] - 1.0;
                    }
                    else
                    {
                        W[n, k] = 1;
                        D[n, k] = 1;
                    }
                }
            }
            var wk = new DenseVector(NK);
            for (int k = 0; k < NK; ++k)
            {
                wk[k] = Math.Sqrt(2.0 * k + 1.0);
            }

            //расчет в один поток
            for (int n = 0; n < NL; ++n)
            {
                ii1 = (int)Math.Ceiling((n + 1.0) / NK) - 1;
                hn = sp(n, NK - 1);
                for (int m = 0; m < NL; ++m)
                {
                    jj1 = (int)Math.Ceiling((m + 1.0) / NK) - 1;
                    hm = sp(m, NK - 1);
                    if (n != m)
                    {
                        Ai[n, m] = int_eigenvaluecorel(ii1, jj1, hn, hm, W, D, X_, zeros_legendre.A0);
                        Ai[n, m] *= wk[hn] * wk[hm];
                        Ae[n, m] = -Ai[n, m];
                    }
                    else
                    {
                        Ai[n, m] = 1;
                        Ae[n, m] = 1;
                    }
                    if ((hn == 0) && (hm == 0)) { Ae[n, m] -= 2.0 * Length_L[jj1]; }
                }
            }

            var Ci = new List<VectorDouble>();
            var Ce = new List<VectorDouble>();
            var F = new DenseVector(NL);
            Matrix<double> ACi = Ai.Inverse();
            Matrix<double> ACe = Ae.Inverse();
            for (int k = 0; k < N; ++k)
            {
                Ci.Add(new DenseVector(NL));
                Ce.Add(new DenseVector(NL));
                for (int n = 0; n < NL; ++n)
                {
                    ii1 = (int)Math.Ceiling((n + 1.0) / NK);
                    hn = sp(n, NK - 1);
                    F[n] = boundary_conditions(ii1, k + 1, hn);
                }
                Ci[k] = ACi * F;
                Ce[k] = -ACe * F;
            }

            return new Coef{
                Ci_ = Ci,
                Ce_ = Ce
            };
        }

        private List<double> lambda_vec(Complex z, int n)
        {
            var zz = legendere.SecondKindVec2(n, z);
            var R = new List<double>();
            for (int i = 0; i < n; ++i)
            {
                R.Add(Math.Sqrt(2 * i + 1.0) * zz[i].Imaginary);
            }
            return R;
        }

        public List<double> barycentric_interier(Complex z)
        {
            int N = P.Count;
            var LL_ = new DenseMatrix(N, NK);
            var zz = new Complex[N];
            double twoPi = 0.5 / Math.PI;
            for (int k = 0; k < N; ++k)
            {
                zz[k] = 2.0 * (z - P[k]) / L[k] - 1.0;
            }
            for (int k = 0; k < N; ++k)
            {
                var ss = lambda_vec(zz[k], NK);
                for (int n = 0; n < NK; ++n)
                {
                    LL_[k, n] = ss[n];
                }
            }
            var zeta = new List<double>(N);
            var ss_ = new DenseVector(N);
            for (int j = 0; j < N; ++j)
            {
                for (int k = 0; k < N; ++k)
                {
                    ss_[k] = Coeff_i[j].Row(k).DotProduct(LL_.Row(k));
                }
                zeta.Add(ss_.Sum() * twoPi);
            }
            return zeta.ToList();
        }

        public List<double> barycentric_exterier(Complex z)
        {
            int N = P.Count;
            var LL_ = new DenseMatrix(N, NK);
            var zz = new Complex[N];
            double twoPi = 0.5 / Math.PI;

            for (int k = 0; k < N; ++k)
            {
                zz[k] = 2.0 * (z - P[k]) / L[k] - 1.0;
            }

            for (int k = 0; k < N; ++k)
            {
                var ss = lambda_vec(zz[k], NK);
                ss[0] += twoPi * Length_L[k];
                for (int n = 0; n < NK; ++n)
                {
                    LL_[k, n] = ss[n];
                }
            }

            var zeta = new List<double>(N);
            var ss_ = new DenseVector(N);

            for (int j = 0; j < N; ++j)
            {
                for (int k = 0; k < N; ++k)
                {
                    ss_[k] = Coeff_e[j].Row(k).DotProduct(LL_.Row(k));
                }
                zeta.Add(ss_.Sum() * twoPi);
            }

            return zeta;
        }

        public bool check_point_in_polygon(Complex p)
        {
            bool result = false;
            int n = P.Count;
            double pImaginary = p.Imaginary;
            double pReal = p.Real;

            for (int i1 = 0, i2 = n - 1; i1 < n; i2 = i1++)
            {
                double P1Imaginary = P[i1].Imaginary;
                double P2Imaginary = P[i2].Imaginary;
                double P1Real = P[i1].Real;
                double P2Real = P[i2].Real;

                if ((P1Imaginary > pImaginary) != (P2Imaginary > pImaginary))
                {
                    double slope = (P2Real - P1Real) / (P2Imaginary - P1Imaginary);
                    double intersectX = slope * (pImaginary - P1Imaginary) + P1Real;
                    if (pReal < intersectX)
                    {
                        result = !result;
                    }
                }
            }
            return result;
        }
    }
}