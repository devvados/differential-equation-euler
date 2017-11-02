using System;
using System.Collections.Generic;
using System.Windows;
using System.Windows.Controls.DataVisualization.Charting;

namespace TempProject
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        //solutions
        public Dictionary<double, double> sourceSolutionsTable = new Dictionary<double, double>(),
                                          derivativeSolutionsTable = new Dictionary<double, double>(),
                                          squaredSolutionsTable = new Dictionary<double, double>(),
                                          squaredDerivativeTable = new Dictionary<double, double>();

        //differential equation parameters
        public static double m, n, t, rhoUpper0, b, g, k, w, H, NSquared, h;
        public static int steps = 0;

        //function to optimize parameters
        double VzDerivativeSquared, VzSquared;

        //differential equation replacement coefficients
        public static double DifEqCoefB, DifEqCoefC;

        public double FunctionToOptimize(double sigma)
        {
            return VzDerivativeSquared * (Math.Pow(sigma, 2) + 4 * (Math.Pow(w, 2))) / (Math.Pow(sigma, 2) * (Math.Pow(n, 2) + Math.Pow(m, 2))) + VzSquared;
        }

        #region defferential equation parameters and their calculation

        /// <summary>
        /// rho0 = Rho0 * exp(-k*z)
        /// </summary>
        /// <param name="z"></param>
        /// <returns></returns>
        public static double Rho0Function(double z)
        {
            return rhoUpper0 * Math.Exp(-k * z);
        }

        /// <summary>
        /// sigma(t) = 2*t
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public static double SigmaFunction(double t)
        {
            return 2 * t;
        }

        /// <summary>
        /// ln(rho0)
        /// </summary>
        /// <returns></returns>
        public static double LnRho0FunctionDerivative()
        {
            return -k;
        }

        /// <summary>
        /// Evaluate initial data
        /// </summary>
        public void EvaluateTransitionalValues()
        {
            NSquared = -g * LnRho0FunctionDerivative();
            h = H / steps;

            DifEqCoefB = LnRho0FunctionDerivative();
            DifEqCoefC = ((Math.Pow(m, 2) + Math.Pow(n, 2)) / (Math.Pow(SigmaFunction(t), 2) - (4 * Math.Pow(w, 2)))) * (Math.Pow(SigmaFunction(t), 2) - NSquared);
        }

        #endregion

        #region differential equations system

        /// <summary>
        /// dy/dz=p
        /// </summary>
        /// <param name="p"></param>
        /// <returns></returns>
        public static double DifferentialEquationPart1(double p)
        {
            return p;
        }

        /// <summary>
        /// dp/dz=C*y-B*p
        /// </summary>
        /// <param name="y"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public static double DifferentialEquationPart2(double y, double p)
        {
            return DifEqCoefC * y - DifEqCoefB * p;
        }

        #endregion

        public MainWindow()
        {
            InitializeComponent();

            sourceSolutionsTable.Add(1, 0);

            TB1Formula.Text = @"V''_z+(\ln\rho_0)'V'_z-\frac{m^2+n^2}{\sigma^2-4w^2}(\sigma^2-N^2)V_z=0";
            TB2Formula.Text = @"V_z|_{z=0}=0, V_z|_{z=-H}=b\sigma";
            TB3Formula.Text = @"E=\int_{-H}^{0}(V'_z)^2dz\cdot\frac{\sigma^2+4w^2}{\sigma^2(n^2+m^2)}+\int_{-H}^{0}V_{z}^{2}dz";
            TB4Formula.Text = @"E=(V'_z)^2\cdot\frac{\sigma^2+4w^2}{\sigma^2(n^2+m^2)}+V_{z}^2";

            //create chart view
            ((LineSeries)LineChart.Series[0]).ItemsSource = sourceSolutionsTable;
        }

        private void BSolveEquation_Click(object sender, RoutedEventArgs e)
        {
            //clear solutions table
            sourceSolutionsTable.Clear();
            derivativeSolutionsTable.Clear();
            squaredSolutionsTable.Clear();
            squaredDerivativeTable.Clear();

            //init values
            InputData();

            //evaluate initial values
            EvaluateTransitionalValues();

            //calculate borders
            double border1 = b * SigmaFunction(t), border2 = 0;

            sourceSolutionsTable = EulerMethod(border1, border2);

            //create chart view
            ((LineSeries)LineChart.Series[0]).ItemsSource = sourceSolutionsTable;

            //fill source function datagrid
            DGSourceFunction.ItemsSource = sourceSolutionsTable;
        }

        private void BSolveIntegral_Click(object sender, RoutedEventArgs e)
        {
            //find derivative
            derivativeSolutionsTable = FindDerivative(sourceSolutionsTable);

            //find squared function
            squaredSolutionsTable = SquareFunction(sourceSolutionsTable);

            //find squared derivative function
            squaredDerivativeTable = SquareFunction(derivativeSolutionsTable);

            //calculate integral
            var resA = SimpsonIntegrate(squaredDerivativeTable, -H, 0, steps);
            var resB = SimpsonIntegrate(squaredSolutionsTable, -H, 0, steps);

            var calculatedE = resA * (Math.Pow(SigmaFunction(t), 2) + 4 * Math.Pow(w, 2)) / ((Math.Pow(SigmaFunction(t), 2)) * (Math.Pow(m, 2) + Math.Pow(n, 2))) + resB;

            //fill function derivative table
            DGFunctionDerivative.ItemsSource = derivativeSolutionsTable;
            //fill derivative squared table
            DGDerivativeFunctionSquared.ItemsSource = squaredDerivativeTable;

            //fill integral value
            TBEValue.Text = Math.Round(calculatedE, 4).ToString();

            //fill source data for optimization
            DGSourceFunctionSquared.ItemsSource = squaredSolutionsTable;
            DGFunctionDerivativeSquared.ItemsSource = squaredDerivativeTable;
        }

        private void BOptimizeFunction_Click(object sender, RoutedEventArgs e)
        {
            var EandSigmaValues = OptimizeFunction();

            //fill datagrid
            DGSigmaAndEValues.ItemsSource = EandSigmaValues;
        }

        /// <summary>
        /// Enter initial data
        /// </summary>
        public void InputData()
        {
            m = Convert.ToDouble(TBM.Text);
            n = Convert.ToDouble(TBN.Text);
            t = Convert.ToDouble(TBT.Text);

            rhoUpper0 = Convert.ToDouble(TBRho.Text);
            b = Convert.ToDouble(TBB.Text);
            w = Convert.ToDouble(TBW.Text);

            g = Convert.ToDouble(TBG.Text);
            k = Convert.ToDouble(TBK.Text);
            H = Convert.ToDouble(TBH.Text);

            steps = Convert.ToInt32(TBSteps.Text);
        }

        /// <summary>
        /// Euler method
        /// </summary>
        /// <param name="b1">left bound</param>
        /// <param name="b2">right bound</param>
        public Dictionary<double, double> EulerMethod(double b1, double b2)
        {
            Dictionary<double, double> temp = new Dictionary<double, double>();
            double K = 0, L = 0, nextB1 = b1, nextB2 = b2;

            for (double i = -H; i <= 0; i += h) 
            {
                nextB1 += h * K;
                nextB2 += h * L;

                K = DifferentialEquationPart1(nextB2);
                L = DifferentialEquationPart2(nextB1, nextB2);

                double tempx = Math.Round(i, 4);
                double tempy = Math.Round(nextB1, 4);

                temp.Add(tempx, tempy);
            }

            return temp;
        }

        /// <summary>
        /// Derivative of table function
        /// </summary>
        /// <param name="d"></param>
        /// <returns></returns>
        public Dictionary<double, double> FindDerivative(Dictionary<double, double> d)
        {
            Dictionary<double, double> temp = new Dictionary<double, double>();
            List<double> xs = new List<double>(d.Keys);
            List<double> ys = new List<double>(d.Values);
            double h = 0.01;
            int n = ys.Count-1;

            var fDerivativeLeftBound = (-3 * ys[0] + 4 * ys[1] - ys[2]) / (2 * h);
            var fDerivativeRightBound = (ys[n - 2] - 4 * ys[n - 1] + 3 * ys[n]) / (2 * h);

            temp.Add(xs[0], fDerivativeLeftBound);

            for(int i = 1; i < xs.Count - 1; i++)
            {
                var fDer = (ys[i + 1] - ys[i - 1]) / (2 * h);
                temp.Add(xs[i], Math.Round(fDer, 4));
            }

            temp.Add(xs[n], fDerivativeRightBound);

            return temp;
        }

        /// <summary>
        /// Squaring table function
        /// </summary>
        /// <param name="d"></param>
        /// <returns></returns>
        public Dictionary<double, double> SquareFunction(Dictionary<double, double> d)
        {
            Dictionary<double, double> temp = new Dictionary<double, double>();

            List<double> xs = new List<double>(d.Keys);
            List<double> ys = new List<double>(d.Values);

            for(int i = 0; i < xs.Count; i++)
            {
                var fSquared = Math.Pow(ys[i], 2);
                temp.Add(xs[i], Math.Round(fSquared, 4));
            }

            return temp;
        }

        public double SimpsonIntegrate(Dictionary<double, double> d, double leftBound, double rightBound, int count)
        {
            List<double> xs = new List<double>(d.Keys);
            List<double> ys = new List<double>(d.Values);
            int n = ys.Count - 1;

            double h = (rightBound - leftBound) / count, result = 0;

            for(int i = 1; i < xs.Count - 2; i++)
            {
                if(i % 2 == 0)
                    result += 4 * ys[i];
                else
                    result += 2 * ys[i];
            }
            result = h / 3 * (result + ys[0] - ys[n]);

            return result; 
        }

        /// <summary>
        /// Newton interpolation
        /// </summary>
        /// <param name="x"></param>
        /// <param name="n"></param>
        /// <param name="MasX"></param>
        /// <param name="MasY"></param>
        /// <param name="step"></param>
        /// <returns></returns>
        public double Newton(double x, int n, List<double> MasX, List<double> MasY, double step)
        {
            double[,] mas = new double[n + 2, n + 1];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < n + 1; j++)
                {
                    if (i == 0)
                        mas[i, j] = MasX[j];
                    else if (i == 1)
                        mas[i, j] = MasY[j];
                }
            }
            int m = n;
            for (int i = 2; i < n + 2; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    mas[i, j] = mas[i - 1, j + 1] - mas[i - 1, j];
                }
                m--;
            }

            double[] dy0 = new double[n + 1];

            for (int i = 0; i < n + 1; i++)
            {
                dy0[i] = mas[i + 1, 0];
            }

            double res = dy0[0];
            double[] xn = new double[n];
            xn[0] = x - mas[0, 0];

            for (int i = 1; i < n; i++)
            {
                double ans = xn[i - 1] * (x - mas[0, i]);
                xn[i] = ans;
                ans = 0;
            }

            int m1 = n + 1;
            int fact = 1;
            for (int i = 1; i < m1; i++)
            {
                fact = fact * i;
                res = res + (dy0[i] * xn[i - 1]) / (fact * Math.Pow(step, i));
            }

            return res;
        }

        /// <summary>
        /// Optimize function -> MIN
        /// </summary>
        /// <returns></returns>
        public Dictionary<double, double> OptimizeFunction()
        {
            Dictionary<double, double> result = new Dictionary<double, double>();
            var keys1 = new List<double>(squaredDerivativeTable.Keys);
            var keys2 = new List<double>(squaredSolutionsTable.Keys);
            double leftBound = -1, rightBound = 0;

            for(int i = 0; i < keys1.Count; i++)
            {                
                VzDerivativeSquared = squaredDerivativeTable[keys1[i]];
                VzSquared = squaredSolutionsTable[keys2[i]];

                var pair = GoldenSectionMethod(leftBound, rightBound);

                if (!result.ContainsKey(pair.Key))
                {
                    result.Add(pair.Key, pair.Value);
                }
                else if (pair.Value < result[pair.Key])
                {
                    result.Remove(pair.Key);
                    result.Add(pair.Key, pair.Value);
                }
            }

            return result;
        }

        /// <summary>
        /// Get E and Sigma values
        /// </summary>
        /// <param name="leftBound"></param>
        /// <param name="rightBound"></param>
        /// <returns></returns>
        public KeyValuePair<double, double> GoldenSectionMethod(double leftBound, double rightBound)
        {
            double tau = (Math.Sqrt(5) - 1) / 2;
            double epsilon = 1e-7;
            double sigma1 = 0, sigma2 = 0, f1 = 0, f2 = 0;
            double a = leftBound, b = rightBound;

            while (b-a > epsilon)
            {
                sigma1 = a + (1 - tau) * (b - a);
                sigma2 = a + tau * (b - a);
                f1 = FunctionToOptimize(sigma1);
                f2 = FunctionToOptimize(sigma2);

                if (f1 > f2)
                {
                    a = sigma1;
                }
                else
                    b = sigma2;
            }
            if((b-a) <= epsilon)
            {
                sigma1 = (b + a) / 2;
                f1 = FunctionToOptimize(sigma1);
            }

            return new KeyValuePair<double, double>(sigma1, f1);
        }
    }
}
