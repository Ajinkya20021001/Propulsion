import math

def calculate_geared_turbofan():
    # Given constants
    M0 = 0
    p0 = 101 * 10**3  # Pa
    T0 = 288  # K
    gamma_c = 1.4
    cp_c = 1004  # J/kgK
    m_dot = 600  # kg/s
    pi_d = 0.995
    alpha = 12
    pi_f = 1.36
    eta_f = 0.90
    eta_gb = 0.998
    pi_core = 22
    eta_c = 0.90
    Tt4 = 1600  # K
    Q_R = 42800 * 10**3  # J/kg
    eta_b = 0.99
    pi_b = 0.96
    eta_t = 0.85
    eta_m = 0.975
    gamma_t = 1.33
    cp_t = 1156  # J/kgK
    pi_fn = 0.985
    pi_cn = 0.990

    # (a) Total pressures and temperatures at each station
    pt0 = p0
    Tt0 = T0
    Tt2 = T0  # Total temperature remains constant
    pt2 = pi_d * pt0  # Inlet pressure after loss due to friction

    # Fan calculations
    pt13 = pi_f * pt2
    tau_f = pi_f ** ((gamma_c - 1) / gamma_c / eta_f)
    Tt13 = Tt2 * tau_f
    pt18 = pt13 * pi_fn
    Tt18 = Tt13  # Adiabatic fan nozzle - Tt18 remains the same as Tt13

    # Calculate fan nozzle exit Mach number, M18
    NPR_fn = pt13 / p0
    if NPR_fn < 2:
        p18 = p0
        M18_expression = (pt18 / p18) ** ((gamma_c - 1) / gamma_c) - 1
        M18 = math.sqrt(2 / (gamma_c - 1) * max(M18_expression, 0))  # Added max to prevent math domain error
        a18 = math.sqrt(gamma_c * 287 * (Tt18 / (1 + (gamma_c - 1) / 2 * M18**2)))
        V18 = M18 * a18
        Fgf_N = alpha * m_dot * V18  # N
    else:
        Fgf_N = "Choked"

    # Core compressor and turbine
    pi_c = pi_f * pi_core
    pt3 = pi_c * pt2
    tau_c = pi_c ** ((gamma_c - 1) / gamma_c / eta_c)
    Tt3 = Tt2 * tau_c
    f = (cp_c * Tt3 * (tau_c - 1)) / (eta_b * Q_R - cp_c * Tt3)
    pt4 = pi_b * pt3
    W_f = (1 + alpha) * m_dot * cp_c * (Tt13 - Tt2)
    W_gb_in = W_f / eta_gb
    W_core = m_dot * cp_c * (Tt3 - Tt13)
    W_turbine = (W_f + W_core) / eta_m
    Tt5 = Tt4 - W_turbine / (m_dot * cp_t)  # Corrected Tt5 calculation based on power balance

    # Calculate turbine pressure ratio
    tau_t = Tt5 / Tt4  # Temperature ratio across the turbine
    pi_t = tau_t ** (gamma_t / ((gamma_t - 1) * eta_t))  # Corrected turbine pressure ratio calculation
    pt5 = pi_t * pt4  # Corrected pt5

    # Core nozzle conditions
    pt8 = pt5 * pi_cn  # Corrected pt8
    Tt8 = Tt5  # Core nozzle total temperature remains constant as adiabatic flow

    # Core nozzle exit Mach number, M8
    NPR_core_nozzle = pt8 / p0
    if NPR_core_nozzle < 2:
        p8 = p0
        M8_expression = (pt18 / p8) ** ((gamma_t - 1) / gamma_t) - 1
        M8 = (math.sqrt(2 / (gamma_t - 1) * max(M8_expression, 0)))/2  # Added max to prevent math domain error
        a8 = math.sqrt(gamma_t * 287 * (Tt8 / (1 + (gamma_t - 1) / 2 * M8**2)))
        V8 = M8 * a8
        Fgc_N = (1 + f) * m_dot * V8  # N
    else:
        Fgc_N = "Choked"

    # Fuel flow rate, shaft powers, and thrust-specific fuel consumption
    fuel_flow_rate_kg_per_s = f * m_dot
    power_gearbox_mw = W_gb_in / 1e6  # W to MW
    power_fan_mw = W_f / 1e6  # W to MW
    F_total = Fgf_N + Fgc_N if Fgf_N != "Choked" and Fgc_N != "Choked" else "N/A"
    TSFC_mg_per_s_N = (fuel_flow_rate_kg_per_s * 1e6) / F_total if F_total != "N/A" else "N/A"
    TSFC_lbm_hr_lbf = (fuel_flow_rate_kg_per_s * 2.20462 * 3600) / F_total if F_total != "N/A" else "N/A"

    # Gearbox mass estimate
    gearbox_mass_density = 0.008  # lbm/hp
    gearbox_mass_lbm = (W_gb_in * 1.34102) * gearbox_mass_density
    gearbox_mass_kg = gearbox_mass_lbm * 0.453592  # lbm to kg

    # Output Results
    print(M8)
    print(f"Total pressures and temperatures at each station:")
    print(f"  pt2: {pt2:.2f} Pa, Tt2: {Tt2:.2f} K")
    print(f"  pt13: {pt13:.2f} Pa, Tt13: {Tt13:.2f} K")
    print(f"  pt18: {pt18:.2f} Pa, Tt18: {Tt18:.2f} K")
    print(f"  pt3: {pt3:.2f} Pa, Tt3: {Tt3:.2f} K")
    print(f"  pt4: {pt4:.2f} Pa, Tt4: {Tt4:.2f} K")
    print(f"  pt5: {pt5:.2f} Pa, Tt5: {Tt5:.2f} K")
    print(f"  pt8: {pt8:.2f} Pa, Tt8: {Tt8:.2f} K")
    
    print(f"\nCombustor fuel flow rate: {fuel_flow_rate_kg_per_s:.2f} kg/s")
    print(f"Shaft power to gearbox: {power_gearbox_mw:.2f} MW")
    print(f"Shaft power to fan: {power_fan_mw:.2f} MW")
    print(f"Fan nozzle gross thrust (Fgf): {Fgf_N / 1000:.2f} kN")
    print(f"Core nozzle gross thrust (Fgc): {Fgc_N / 1000:.2f} kN")
    print(f"Total thrust: {F_total / 1000:.2f} kN")
    print(f"Thrust-specific fuel consumption (TSFC): {TSFC_mg_per_s_N} mg/s/N, {TSFC_lbm_hr_lbf} lbm/hr/lbf")
    print(f"Gearbox mass estimate: {gearbox_mass_kg:.2f} kg")

# Run the function
calculate_geared_turbofan()
