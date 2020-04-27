/* Preciso calcular entalpia, entropia e fazer uma função que receba composição, 
temperatura e pressão e retorne as fases em equilibrio e as composições das fases.
*/
using System;

namespace BayerProcess
{
    class functions
    {
        public double BPE_Dewey( double Molality, double Pressure_kPa )
        {
            double conv = 1;
            double BPE = 0;
            double temp_final;
            double water_BP = 100; // This water boiling point needs to be calculated by some water property package based on the pressure given
            double loops = 0;

            while (conv > 1e-8 && loops < 1000)
            {
                temp_final = water_BP + BPE;
                BPE = 1.82e-3 + 0.55379 * Math.Pow( Molality / 10 , 7 );
                BPE += 4.0625e-3 * Molality * temp_final;
                BPE += 1 / temp_final * ( -286.66 * Molality + 29.919 * Math.Pow( Molality ,2 ) +0.6228 * Math.Pow( Molality , 3 ) );
                BPE -= 32.647e-3 * Molality * Math.Pow( Molality * temp_final * 1e-3 , 2 ) ;
                BPE += ( Math.Pow( temp_final * 1e-3 , 5 ) ) * ( 5.9705 * Molality - ( 0.57532 * Math.Pow( Molality ,2 ) ) + ( 0.10417 * Math.Pow( Molality, 3 ) ) );

                conv = temp_final - ( water_BP + BPE );
                loops++;
            }

            return BPE;
        }

        public double BPE_Adamson( double total_soda_Na2Ogpl, double Pressure_kPa)
        {
            double conv = 1;
            double BPE = 0;
            double temp_final;
            double water_BP = 100; // This water boiling point needs to be calculated by some water property package based on the pressure given
            int loops = 0;

            while (conv > 1e-8 && loops < 1000 )
            {
                temp_final = water_BP + BPE;
                BPE = 7.642857e-3 + 6.184282e-3 * total_soda_Na2Ogpl + 2.92857e-5 * temp_final; // independent and linear terms
                BPE += 1.0957e-4 * Math.Pow( total_soda_Na2Ogpl, 2 ) - 3.80952e-8 * Math.Pow( temp_final, 2 ) + 2.08801e-4 * total_soda_Na2Ogpl * temp_final; // quadratic terms
                BPE += -8.61985e-10 * Math.Pow( total_soda_Na2Ogpl, 3 ) - 8.61985e-10 * Math.Pow( temp_final, 3) + 1.7316e-10 * total_soda_Na2Ogpl * Math.Pow( temp_final, 2 )  - 2.49763e-7 * Math.Pow( total_soda_Na2Ogpl, 2) *temp_final;

                conv = temp_final - ( water_BP + BPE );
                loops++;
            }   

            return BPE;
        }

        public double Cp_Langa ( double total_caustic_Na2CO3gpl, double alumina_Al2O3gpl, double temperature )
        {
            /* 
            Function programmed by Daniel Rodrigues on 24/04/2020 for an alumina DWSIM package based on Langa's (1985) Cp correlation  
            */

            double Cp, CP1, CP2;

            CP1 = 0.99639 - 3.90998e-4 * total_caustic_Na2CO3gpl - 5.3832e-4 * alumina_Al2O3gpl + 2.46493e-7 * Math.Pow( total_caustic_Na2CO3gpl, 2 ) + 5.7186e-7 * total_caustic_Na2CO3gpl * alumina_Al2O3gpl;
            CP2 = -1.51278e-4 - 1.86581e-7 * alumina_Al2O3gpl - 1.07766e-7 * total_caustic_Na2CO3gpl;
            Cp = 4.184 * ( CP1 + temperature * ( CP2 + 2.1464e-6 * temperature ) );

            return Cp; // kJ/kg
        }

        public double Entropy_liquor ( double total_caustic_Na2CO3gpl, double alumina_Al2O3gpl, double temperature )
        {
            /* 
            Function programmed by Daniel Rodrigues on 27/04/2020 for an alumina DWSIM package based on Langa's (1985) Cp correlation.
            S(T) = ref_ent + integral (Cp/Temp , Temp)  = 4.184 * ( ln(T/Tref) + K2*(T-Tref) + K3*(T^2 - Tref^2)/2 )
            assuming solution entropy at 25C = 0 
            */

            double Entropy, CP2;

            CP2 = -1.51278e-4 - 1.86581e-7 * alumina_Al2O3gpl - 1.07766e-7 * total_caustic_Na2CO3gpl;
            Entropy = 4.184 * ( Math.Log(temperature / 25) + CP2 * ( temperature - 25 ) + 2.1464e-6 /2 * ( Math.Pow( temperature, 2 ) - Math.Pow( 25, 2 ) ) );

            return Entropy; //kJ/kg C
        }

        public double Enthalpy_liquor ( double total_caustic_Na2CO3gpl, double alumina_Al2O3gpl, double temperature )
        {   
            /* 
            Function programmed by Daniel Rodrigues on 27/04/2020 for an alumina DWSIM package based on Langa's (1985) Cp correlation.
            H(T) = ref_ent + integral (Cp, Temp)  
            assuming solution enthalpy at 25C = 0 
            */

            double Enthalpy, CP1, CP2;

            CP1 = 0.99639 - 3.90998e-4 * total_caustic_Na2CO3gpl - 5.3832e-4 * alumina_Al2O3gpl + 2.46493e-7 * Math.Pow( total_caustic_Na2CO3gpl, 2 ) + 5.7186e-7 * total_caustic_Na2CO3gpl * alumina_Al2O3gpl;
            CP2 = -1.51278e-4 - 1.86581e-7 * alumina_Al2O3gpl - 1.07766e-7 * total_caustic_Na2CO3gpl;
            Enthalpy = 4.184 * ( CP1 * temperature  + CP2  / 2 * Math.Pow( temperature, 2 ) + 2.1464e-6 / 3 * Math.Pow( temperature, 3 ) );

            return Enthalpy; //kJ/kg
        } 

        public double Liquor_density_mass_fraction ( double total_soda_Na2O_mf, double alumina_Al2O3_mf, double temperature)
        {
            /*
            Function programmed by Daniel Rodrigues on 27/04/2020 for an alumina DWSIM package based on Mulloy-Donaldson's density correlation.
            */
            double LSG_T, LSG_25;

            LSG_25 = 0.982;
            LSG_25 +=  1.349855e-2 * total_soda_Na2O_mf - 2.4948e-4 * Math.Pow( total_soda_Na2O_mf, 2 ) + 2.73e-6 * Math.Pow( total_soda_Na2O_mf, 3);
            LSG_25 += 2.08035e-3 * alumina_Al2O3_mf - 7.28e-6 * Math.Pow( alumina_Al2O3_mf, 2 );
            LSG_25 += 3.3367e-4 * total_soda_Na2O_mf * alumina_Al2O3_mf;

            LSG_T = LSG_25 * ( 1 - ( 5.021858e-4 * 0.85 * ( temperature - 25 ) ) - ( 1.1881e-6 * 0.85 * Math.Pow(  temperature - 25, 2 ) ) );

            return LSG_T;
        }

        public double Liquor_density_conc ( double total_soda_Na2Ogpl, double alumina_Al2O3gpl, double temperature )
        {
            /*
            Function programmed by Daniel Rodrigues on 27/04/2020 for an alumina DWSIM package based on Mulloy-Donaldson's density correlation.
            */
            double LSG_T, LSG_25, initDensity, total_soda_Na2O_mf = 0, alumina_Al2O3_mf = 0;
            int counter = 0;

            LSG_25 = 1.25;
            initDensity = 1.0;
	        counter = 0;
	        while ( Math.Abs( LSG_25 - initDensity ) > 1e-8 && counter < 1000 ) 
            {
                initDensity = LSG_25;
	            total_soda_Na2O_mf = total_soda_Na2Ogpl / initDensity / 10;
	            alumina_Al2O3_mf = alumina_Al2O3gpl / initDensity / 10;
	            LSG_25 = Liquor_density_mass_fraction ( total_soda_Na2O_mf, alumina_Al2O3_mf, 25 );
    
	            counter++;
	        
            }   
	    
	        LSG_T = Liquor_density_mass_fraction ( total_soda_Na2O_mf, alumina_Al2O3_mf, temperature );
	
        return LSG_T;
        }

        public double[] flash_calc ( double liquor_in_temp, double vessel_P_kPa , double liquor_in_Al2O3gpl_25C, double liquor_in_TC_Na2CO3gpl_25C, double liquor_in_total_soda_Na2CO3gpl_25C )
        {
            /*
            Function programmed by Daniel Rodrigues on 27/04/2020 for an alumina DWSIM package.
            */
            double new_evaporated_mass_fraction , liquor_in_density, liquor_out_density, liquor_in_total_soda_Na2CO3_mf, liquor_in_TC_Na2CO3_mf, liquor_in_Al2O3_mf;
            double liquor_out_Al2O3_mf, liquor_out_total_soda_Na2O_mf, liquor_out_TC_Na2CO3_mf, liquor_out_Al2O3gpl = 0 , liquor_out_TC_Na2CO3_gpl = 0, liquor_out_total_soda_Na2CO3gpl = 0;
            double enth_in, enth_out;
            double conv  = 1, final_boiling_point = 100, evaporated_mass_fraction = 0;
            double water_BP = 100, superheated_steam_enthalpy = 1000; // this is probably going to come in somehow from DWSIM
            double[] result = new double[4];
            int counter = 0;

            liquor_in_density = Liquor_density_conc( liquor_in_total_soda_Na2CO3gpl_25C, liquor_in_Al2O3gpl_25C, 25 );
            liquor_in_total_soda_Na2CO3_mf = liquor_in_total_soda_Na2CO3gpl_25C /  liquor_in_density / 1000;
            liquor_in_TC_Na2CO3_mf = liquor_in_TC_Na2CO3gpl_25C / liquor_in_density / 1000;
            liquor_in_Al2O3_mf = liquor_in_Al2O3gpl_25C / liquor_in_density / 1000;

            while ( conv > 1e-8 && counter < 1000 )
            {
                liquor_out_Al2O3_mf = liquor_in_Al2O3_mf / ( 1 - evaporated_mass_fraction );
                liquor_out_total_soda_Na2O_mf = liquor_in_total_soda_Na2CO3_mf / ( 1 - evaporated_mass_fraction );
                liquor_out_TC_Na2CO3_mf = liquor_in_TC_Na2CO3_mf / ( 1 - evaporated_mass_fraction );

                liquor_out_density = Liquor_density_mass_fraction( liquor_out_total_soda_Na2O_mf * 100, liquor_out_Al2O3_mf * 100, 25 );
                liquor_out_Al2O3gpl = liquor_out_Al2O3_mf * liquor_out_density * 1000;
                liquor_out_TC_Na2CO3_gpl = liquor_out_TC_Na2CO3_mf * liquor_out_density * 1000;
                liquor_out_total_soda_Na2CO3gpl = liquor_out_total_soda_Na2O_mf * liquor_out_density * 1000; 

                final_boiling_point = water_BP + BPE_Adamson( liquor_out_total_soda_Na2CO3gpl, vessel_P_kPa );
                enth_in = Enthalpy_liquor( liquor_in_Al2O3gpl_25C, liquor_in_Al2O3gpl_25C, liquor_in_temp );
                enth_out = Enthalpy_liquor( liquor_out_TC_Na2CO3_gpl, liquor_out_Al2O3gpl, final_boiling_point );
                new_evaporated_mass_fraction = ( enth_in - enth_out ) / ( superheated_steam_enthalpy - enth_out );
        
                conv = Math.Abs( evaporated_mass_fraction - new_evaporated_mass_fraction );
                evaporated_mass_fraction = new_evaporated_mass_fraction;
                counter++;
            }

            result[0] = final_boiling_point;
            result[1] = liquor_out_Al2O3gpl;
            result[2] = liquor_out_TC_Na2CO3_gpl;
            result[3] = liquor_out_total_soda_Na2CO3gpl;
    
            return result;
        }
    }
}
