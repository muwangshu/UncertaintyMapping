using System;
using System.Collections.Generic;
using System.IO;


namespace SegmentGeneration
{
    struct Brk
    {
        public double value;
        public int idx;
    }
    class Program 
    {
        static void Main(string[] args)
        {
            args = new string[6] { @"..\..\..\..\Data\MHI_2014_us_county.csv", "outtest.txt", "5", "0", "80", "80" };
            string infile = args[0];
            string outfile = args[1];
            int p = Convert.ToInt32(args[2]);
            int interval = Convert.ToInt32(args[3]);
            int rate = Convert.ToInt32(args[4]);
            int confi = Convert.ToInt32(args[5]);
            Normal(infile, outfile,p, ((double)rate) / 100, ((double)confi / 100), interval);
        }

        static Brk[] GenerateBreaks(double[] estimates)
        {
            double[] tbrks = new double[estimates.Length - 1];
            int j = 0,i=0;
            while(j<estimates.Length-1)
            {
                if (estimates[j] < estimates[j + 1])
                {
                    tbrks[i] = (double)(estimates[j] + estimates[j + 1]) / 2;
                    i++;
                    j++;
                }
                else
                {
                    j++;
                }
            }
            Brk[] brks = new Brk[i + 2];
            int n = i;
            for (i = 0; i < n; i++)
                brks[i + 1].value = tbrks[i];

            brks[0].value = double.MinValue;
            brks[0].idx = -1;
            brks[n + 1].value = double.MaxValue;
            brks[n + 1].idx = estimates.Length - 1;

            int ii = 0;
            for (int t = 1; t <= n; t++)
            {
                while (estimates[ii] <= brks[t].value) ii++;

                ii--;
                brks[t].idx = ii;
            }

            return brks;
        }
        static Brk[] GenerateBreaks(double[] estimates, int interval)
        {
            double start = estimates[0];
            double end = estimates[estimates.Length - 1];
            double s = Math.Ceiling(start / interval) * interval;
            double e = Math.Floor(end / interval) * interval;
            int n = (int)(e - s) / interval + 1;
            if (e == end) n = n - 1;
            Brk[] brks = new Brk[n+2];
            for (int i=1;i<=n;i++)
            {
                brks[i].value = s + interval * (i-1);
            }
            brks[0].value = double.MinValue;
            brks[0].idx = -1;
            brks[n + 1].value = double.MaxValue;
            brks[n + 1].idx = estimates.Length-1;

            int ii = 0;
            for (int t = 1; t<=n;t++)
            {
                while (estimates[ii] <=brks[t].value) ii++;

                ii--;
                brks[t].idx = ii;
            }

            return brks;
        }

        static void Normal(string infile, string outfile,int p,double rate,double confi,int interval)
        {

            var watch = System.Diagnostics.Stopwatch.StartNew();

            
            string[] content = File.ReadAllLines(infile);
            int count = content.Length;

            double[] estimate = new double[count];
            double[] error = new double[count];
            for (int i = 0; i < count; i++)
            {
                string[] dataline = content[i].Split(',');
                estimate[i] = Convert.ToDouble(dataline[0]);
                error[i] = Math.Pow(Convert.ToDouble(dataline[1]) / 1.645, 2);
            }


            Brk[] brks;
            if (interval == 0)
                brks = GenerateBreaks(estimate);
            else
                brks = GenerateBreaks(estimate, interval);

            int bkcount = brks.Length;
            double[,] CPD = CalculateCumulativePD(estimate, error, brks);
            double[,] ECPD = CalculateCumulativePD(estimate, error);
            double[,] objvals = new double[bkcount, bkcount];
            for (int i = 0; i < bkcount; i++)
                for (int j = 0; j < bkcount; j++)
                    objvals[i, j] = double.MaxValue;
            double[,] optimals = new double[bkcount, bkcount];

            for (int i = 0; i < bkcount - 1; i++)
                for (int j = i + 1; j < bkcount; j++)
                {
                    //if (Validate(estimate, error, brks[i].idx + 1, brks[j].idx, brks[i].value, brks[j].value, rate, confi) == true)
                    if ((interval!=0 && Validate(CPD, brks[i].idx + 1, brks[j].idx,i,j,rate,confi)==true) ||
                    ( interval==0 && ValidateExpand(CPD,ECPD,estimate,brks, brks[i].idx + 1, brks[j].idx, i, j, rate, confi)) == true)
                    {
                        double upper = 0, down = 0;
                        for (int k = brks[i].idx + 1; k <= brks[j].idx; k++)
                        {
                            upper += estimate[k] / error[k];
                            down += 1 / error[k];
                        }
                        double optimal = upper / down;
                        optimals[i, j] = optimal;
                        double objval = 0;
                        for (int k = brks[i].idx + 1; k <= brks[j].idx; k++)
                        {
                            objval += (optimal - estimate[k]) * (optimal - estimate[k]) / error[k];
                        }
                        objvals[i, j] = objval;
                    }
                    else
                    {
                        objvals[i, j] = double.MaxValue;
                    }
                }


            double[,] dpspace = new double[p, bkcount];
            int[,] track = new int[p, bkcount];
            for (int i = 0; i < bkcount; i++)
                dpspace[0, i] = objvals[0, i];
            for (int i = 1; i < p; i++)
            {
                for (int j = i+1; j < bkcount; j++)
                {
                    double min = double.MaxValue;
                    int mink = -1;
                    for (int k = 0; k < j; k++)
                    {
                        double val = dpspace[i - 1, k] + objvals[k, j];
                        if (val < min) { min = val; mink = k; }
                    }
                    dpspace[i, j] = min;
                    track[i, j] = mink;
                }
            }

            string filename = outfile;
            if (File.Exists(filename)) File.Delete(filename);
            int jj = bkcount - 1;
            for (int i = p - 1; i >= 0; i--)
            {
                int k = jj;
                jj = track[i, jj];
                double optimal = optimals[jj, k];
                double objval = objvals[jj, k];
                File.AppendAllText(filename, string.Format("{0},{1},{2},{3},{4}\r\n", jj, brks[jj].idx,brks[jj].value, optimal, objval));
            }
            File.AppendAllText(filename, string.Format("{0}\r\n", dpspace[p - 1, bkcount - 1]));
            watch.Stop();
            var elapsedSec = watch.ElapsedMilliseconds/1000;
            File.AppendAllText(filename, string.Format("{0}\r\n", elapsedSec));
        }
        static double[,] CalculateCumulativePD(double[] estimate, double[] error,Brk[] brks)
        {
            double[,] CPD = new double[estimate.Length, brks.Length];

            for (int i=0;i<estimate.Length;i++)
                for (int j=0;j<brks.Length;j++)
                {
                    CPD[i,j]= NormalDistributionConfidenceCalculator.CumDensity((brks[j].value - estimate[i]) / Math.Sqrt(error[i]));
                }
            return CPD;
        }
        static double[,] CalculateCumulativePD(double[] estimate, double[] error)
        {
            double[,] CPD = new double[estimate.Length, estimate.Length];

            for (int i = 0; i < estimate.Length; i++)
                for (int j = 0; j < estimate.Length; j++)
                {
                    CPD[i, j] = NormalDistributionConfidenceCalculator.CumDensity((estimate[j] - estimate[i]) / Math.Sqrt(error[i]));
                }
            return CPD;
        }
        static bool Validate(double[,] CPD,int ui, int uj,int brki,int brkj, double rate, double confidencelevel)
        {
            int count = CPD.GetLength(0);
            double[] leftprob = new double[count];
            double[] rightprob = new double[count];

            if (uj == int.MaxValue) uj = count - 1;
            if (ui > uj) return false;
            if (ui > 0)
            {
                for (int i = ui; i <= uj; i++)
                {
                    leftprob[i] = CPD[i, brki];
                }
            }
            if (uj == count - 1)
            {
                for (int i = ui; i <= uj; i++)
                {
                    rightprob[i] = 1;
                }


            }
            else
            {
                for (int i = ui; i <= uj; i++)
                {
                    rightprob[i] = CPD[i, brkj];
                }
            }
            int valid = 0;
            for (int i = ui; i <= uj; i++)
            {
                if (rightprob[i] - leftprob[i] >= confidencelevel) valid++;
            }

            if (valid >= (uj - ui + 1) * rate) return true;
            else return false;
        }
       
        static bool ValidateExpand(double[,] CPD, double[,] ECPD, double[] estimate, Brk[] brks, int ui, int uj, int brki, int brkj, double rate, double confidencelevel)
        {
            // expanded range: left: max(brks[brki-1].value,estimate[ui-1])  right: min(brks[brkj+1].value,estimate[j+1])

            int count = CPD.GetLength(0);
            double[] leftprob = new double[count];
            double[] rightprob = new double[count];
            double[,] t;
            int ti,tj;
            if (uj == int.MaxValue) uj = count - 1;
            if (ui > uj) return false;
            if (ui > 0)
            {
                if (brki==0)
                {
                    t = CPD;
                    ti = brki;
                }
                else if (brks[brki-1].value>estimate[ui-1])
                {
                    t = CPD;
                    ti = brki-1;
                }
                else
                {
                    t = ECPD;
                    ti = ui - 1;
                }
                for (int i = ui; i <= uj; i++)
                {
                    leftprob[i] = t[i, ti];
                }
            }

            if (uj == count - 1)
            {
                for (int i = ui; i <= uj; i++)
                {
                    rightprob[i] = 1;
                }
            }
            else
            {
                if (brkj==brks.Length-1)
                {
                    t = CPD;
                    tj = brkj;
                }
                else if (brks[brkj+ 1].value < estimate[uj + 1])
                {
                    t = CPD;
                    tj = brkj+1;
                }
                else
                {
                    t = ECPD;
                    tj = uj + 1;
                }
                for (int i = ui; i <= uj; i++)
                {
                    rightprob[i] = t[i, tj];
                }
            }

            int valid = 0;
            for (int i = ui; i <= uj; i++)
            {
                if (rightprob[i] - leftprob[i] >= confidencelevel) valid++;
            }

            if (valid >= (uj - ui + 1) * rate) return true;
            else return false;
        }
    }
}
