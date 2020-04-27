using System;

namespace BayerProcess
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
            double[] results;
            double density_gpl, density_mf;
            functions fs = new functions();
            
            density_gpl = fs.Liquor_density_conc( 310, 150, 25 );
            density_mf = fs.Liquor_density_mass_fraction ( 31 / density_gpl, 15/density_gpl, 25 );
            results = fs.flash_calc( 150, 101.125, 150, 300, 310 );
            Console.WriteLine(results[0]); // final temperature of superheated steam and liquor leaving the flash
            Console.WriteLine(results[1]); // new Al2O3 conc in the liquor leaving the flash
            Console.WriteLine(results[2]); // new TC conc in the liquor leaving the flash in Na2CO3 gpl
            Console.WriteLine(results[3]); // new TS conc in the liquor leaving the flash in Na2CO3 gpl

        }
    }
}
