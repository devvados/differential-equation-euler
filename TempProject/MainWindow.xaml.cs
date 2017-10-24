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
        public Dictionary<double, double> sourceSolutionsTable = new Dictionary<double, double>();

        //differential equation parameters
        public static double m, n, t, rhoUpper0, b, g, k, w, H, NSquared, h;
        public static int steps = 0;

        //differential equation replacement coefficients
        public static double DifEqCoefB, DifEqCoefC;

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

            TB1Formula.Text = @"V''_z+(\ln\rho_0)'V'_z-\frac{m^2+n^2}{\sigma^2-4w^2}(\sigma^2-N^2)V_z=0";
            TB2Formula.Text = @"V_z|_{z=0}=0, V_z|_{z=-H}=b\sigma";
            TB3Formula.Text = @"E=\int_{-H}^{0}(V'_z)^2dz\cdot\frac{\sigma^2+4w^2}{\sigma^2(n^2+m^2)}+\int_{-H}^{0}V_{z}^{2}dz";
        }

        private void BSolveEquation_Click(object sender, RoutedEventArgs e)
        {
            //clear solutions table
            sourceSolutionsTable.Clear();

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
            var derivativeSolutionsTable = FindDerivative(sourceSolutionsTable);

            //find squared function
            var squaredSolutionsTable = SquareFunction(sourceSolutionsTable);

            //find squared derivative function
            var squaredDerivativeTable = SquareFunction(derivativeSolutionsTable);

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
    }
}
