import math

class TurbofanEngine:
    def __init__(self):
        # Constants
        self.R = 287  # Gas constant for air, J/kg.K

    def calculate(self, inputs):
        """
        Calculate all parameters for the geared turbofan engine.
        """
        # Initialize input parameters
        self.m0 = inputs['m0']  # Mass flow rate (kg/s)
        self.p0 = inputs['p0']  # Ambient pressure (kPa)
        self.T0 = inputs['T0']  # Ambient temperature (K)
        self.gamma_c = inputs['gamma_c']  # Cold section specific heat ratio
        self.gamma_t = inputs['gamma_t']  # Hot section specific heat ratio
        self.cpc = inputs['cpc']  # Cold section specific heat (J/kg.K)
        self.cpt = inputs['cpt']  # Hot section specific heat (J/kg.K)
        self.pi_d = inputs['pi_d']  # Inlet pressure recovery
        self.alpha = inputs['alpha']  # Bypass ratio
        self.pi_f = inputs['pi_f']  # Fan pressure ratio
        self.ef = inputs['ef']  # Fan polytropic efficiency
        self.eta_gb = inputs['eta_gb']  # Gearbox efficiency
        self.pi_c = inputs['pi_c']  # Core compressor pressure ratio
        self.ec = inputs['ec']  # Compressor polytropic efficiency
        self.T_t4 = inputs['T_t4']  # Turbine inlet temperature (K)
        self.QR = inputs['QR']  # Fuel heating value (kJ/kg)
        self.eta_b = inputs['eta_b']  # Combustor efficiency
        self.pi_b = inputs['pi_b']  # Combustor pressure ratio
        self.et = inputs['et']  # Turbine polytropic efficiency
        self.eta_m = inputs['eta_m']  # Mechanical efficiency
        self.pi_fn = inputs['pi_fn']  # Fan nozzle pressure ratio
        self.pi_cn = inputs['pi_cn']  # Core nozzle pressure ratio

        results = {}
        
        # Station 0 (Ambient)
        results['p0'] = self.p0
        results['T0'] = self.T0
        results['pt0'] = self.p0  # Total pressure at ambient
        results['Tt0'] = self.T0  # Total temperature at ambient
        
        # Station 2 (After inlet)
        results['pt2'] = results['pt0'] * self.pi_d
        results['Tt2'] = results['Tt0']
        
        # Station 13 (Fan exit - bypass)
        results['pt13'] = results['pt2'] * self.pi_f
        results['Tt13'] = results['Tt2'] * (self.pi_f ** ((self.gamma_c - 1)/(self.gamma_c * self.ef)))
        
        # Station 3 (After compressor)
        results['pt3'] = results['pt13'] * self.pi_c
        results['Tt3'] = results['Tt13'] * (self.pi_c ** ((self.gamma_c - 1)/(self.gamma_c * self.ec)))
        
        # Station 4 (After combustor)
        results['pt4'] = results['pt3'] * self.pi_b
        results['Tt4'] = self.T_t4
        
        # Calculate fuel flow rate
        m_core = self.m0 / (1 + self.alpha)
        cp_avg = (self.cpc + self.cpt) / 2
        results['fuel_flow'] = (m_core * cp_avg * (self.T_t4 - results['Tt3'])) / (self.QR * 1000 * self.eta_b)
        
        # Calculate turbine temperature ratio
        work_compressor = m_core * self.cpc * (results['Tt3'] - results['Tt13'])
        turbine_work = work_compressor / self.eta_m
        results['Tt5'] = results['Tt4'] - turbine_work / (m_core * self.cpt)
        results['pt5'] = results['pt4'] * (results['Tt5']/results['Tt4']) ** (self.gamma_t/((self.gamma_t-1)*self.et))
        
        # Calculate nozzle parameters
        # Fan nozzle
        pr_crit_cold = (2/(self.gamma_c + 1)) ** (self.gamma_c/(self.gamma_c - 1))
        pr_fn = (results['pt13'] * self.pi_fn) / self.p0
        
        if pr_fn > pr_crit_cold:
            results['M18'] = 1.0
            T18_T13 = 2 / (self.gamma_c + 1)
        else:
            results['M18'] = math.sqrt(2/(self.gamma_c-1) * ((pr_fn**((self.gamma_c-1)/self.gamma_c)) - 1))
            T18_T13 = 1 / (1 + ((self.gamma_c-1)/2) * results['M18']**2)
        
        results['T18'] = results['Tt13'] * T18_T13
        results['V18'] = results['M18'] * math.sqrt(self.gamma_c * self.R * results['T18'])
        
        # Core nozzle
        pr_crit_hot = (2/(self.gamma_t + 1)) ** (self.gamma_t/(self.gamma_t - 1))
        pr_cn = (results['pt5'] * self.pi_cn) / self.p0
        
        if pr_cn > pr_crit_hot:
            results['M8'] = 1.0
            T8_T5 = 2 / (self.gamma_t + 1)
        else:
            results['M8'] = math.sqrt(2/(self.gamma_t-1) * ((pr_cn**((self.gamma_t-1)/self.gamma_t)) - 1))
            T8_T5 = 1 / (1 + ((self.gamma_t-1)/2) * results['M8']**2)
        
        results['T8'] = results['Tt5'] * T8_T5
        results['V8'] = results['M8'] * math.sqrt(self.gamma_t * self.R * results['T8'])
        
        # Calculate thrust
        # Fan thrust
        m_bypass = self.m0 * self.alpha / (1 + self.alpha)
        results['Fg_fan'] = m_bypass * results['V18']
        
        # Core thrust
        m_core_out = m_core + results['fuel_flow']
        results['Fg_core'] = m_core_out * results['V8']
        
        # Calculate TSFC
        total_thrust = results['Fg_fan'] + results['Fg_core']
        results['TSFC'] = (results['fuel_flow'] * 1e6) / total_thrust  # mg/N/s
        
        # Calculate gearbox power
        fan_power = m_bypass * self.cpc * (results['Tt13'] - self.T0)
        results['gearbox_power_in'] = fan_power / self.eta_gb  # Watts
        results['gearbox_power_out'] = fan_power  # Watts
        
        # Calculate gearbox mass
        hp_out = results['gearbox_power_out'] / 745.7  # Convert W to hp
        results['gearbox_mass'] = 0.008 * hp_out * 0.4536  # Convert to kg
        
        return results

def convert_units(results):
    """Convert results to various units for display"""
    conversions = {
        'fuel_flow_lbm_s': results['fuel_flow'] * 2.205,
        'gearbox_power_in_MW': results['gearbox_power_in'] / 1e6,
        'gearbox_power_in_hp': results['gearbox_power_in'] / 745.7,
        'gearbox_power_out_MW': results['gearbox_power_out'] / 1e6,
        'gearbox_power_out_hp': results['gearbox_power_out'] / 745.7,
        'Fg_fan_kN': results['Fg_fan'] / 1000,
        'Fg_fan_lbf': results['Fg_fan'] * 0.225,
        'Fg_core_kN': results['Fg_core'] / 1000,
        'Fg_core_lbf': results['Fg_core'] * 0.225,
        'TSFC_lbm_hr_lbf': results['TSFC'] * 0.0036,  # Convert from mg/N/s to lbm/hr/lbf
        'gearbox_mass_lbm': results['gearbox_mass'] * 2.205
    }
    return conversions

def main():
    # Default inputs based on problem statement
    default_inputs = {
        'm0': 600,  # kg/s
        'p0': 101,  # kPa
        'T0': 288,  # K
        'gamma_c': 1.4,
        'gamma_t': 1.33,
        'cpc': 1004,  # J/kg.K
        'cpt': 1156,  # J/kg.K
        'pi_d': 0.995,
        'alpha': 12,
        'pi_f': 1.36,
        'ef': 0.90,
        'eta_gb': 0.998,
        'pi_c': 22,
        'ec': 0.90,
        'T_t4': 1600,  # K
        'QR': 42800,  # kJ/kg
        'eta_b': 0.99,
        'pi_b': 0.96,
        'et': 0.85,
        'eta_m': 0.975,
        'pi_fn': 0.985,
        'pi_cn': 0.990
    }
    
    # Create engine object
    engine = TurbofanEngine()
    
    # Get user inputs or use defaults
    print("Enter values (press Enter to use default):")
    user_inputs = {}
    for key, default in default_inputs.items():
        while True:
            try:
                value = input(f"{key} [{default}]: ")
                if value == "":
                    user_inputs[key] = default
                else:
                    user_inputs[key] = float(value)
                break
            except ValueError:
                print("Please enter a valid number")
    
    # Calculate results
    results = engine.calculate(user_inputs)
    converted = convert_units(results)
    
    # Display results
    print("\nResults:")
    print(f"Fuel flow rate: {results['fuel_flow']:.2f} kg/s ({converted['fuel_flow_lbm_s']:.2f} lbm/s)")
    print(f"Gearbox input power: {converted['gearbox_power_in_MW']:.2f} MW ({converted['gearbox_power_in_hp']:.2f} hp)")
    print(f"Gearbox output power: {converted['gearbox_power_out_MW']:.2f} MW ({converted['gearbox_power_out_hp']:.2f} hp)")
    print(f"Fan nozzle exit Mach number: {results['M18']:.3f}")
    print(f"Core nozzle exit Mach number: {results['M8']:.3f}")
    print(f"Fan nozzle gross thrust: {converted['Fg_fan_kN']:.2f} kN ({converted['Fg_fan_lbf']:.2f} lbf)")
    print(f"Core nozzle gross thrust: {converted['Fg_core_kN']:.2f} kN ({converted['Fg_core_lbf']:.2f} lbf)")
    print(f"TSFC: {results['TSFC']:.2f} mg/N/s ({converted['TSFC_lbm_hr_lbf']:.3f} lbm/hr/lbf)")
    print(f"Estimated gearbox mass: {results['gearbox_mass']:.2f} kg ({converted['gearbox_mass_lbm']:.2f} lbm)")
    
    # Print additional temperatures and pressures
    print("\nStation Properties:")
    for station in [2, 13, 3, 4, 5]:
        print(f"\nStation {station}:")
        if f'pt{station}' in results:
            print(f"Total Pressure: {results[f'pt{station}']:.2f} kPa")
        if f'Tt{station}' in results:
            print(f"Total Temperature: {results[f'Tt{station}']:.2f} K")

if __name__ == "__main__":
    main()

    