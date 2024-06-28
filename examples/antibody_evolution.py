import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from emodlib.malaria import IntrahostComponent
# IntrahostComponent.set_params()

def age_based_surface_area(age_in_years):
    # piecewise linear rising from birth to age 2
    # and then shallower slope to age 20
    newborn_risk = 0.07
    two_year_old_risk = 0.23
    if age_in_years < 2:
        return newborn_risk + age_in_years * (two_year_old_risk - newborn_risk) / 2

    if age_in_years < 20:
        return two_year_old_risk + (age_in_years - 2) * (1 - two_year_old_risk) / ((20 - 2))

    return 1.
def predict_emod_pfemp1_variant_fraction(age_in_years, relative_biting_rate, daily_sim_eir):
    # Predict equilibrium immunity level that EMOD would give based on age, relative biting rate, and daily simulated EIR
    # From logistic function fit across EMOD sweep
    # see C:\Users\joshsu\OneDrive - Bill & Melinda Gates Foundation\Code\emod-network-testing\analysis\240520\summarize_emod_infections.ipynb

    individual_daily_eir = daily_sim_eir * relative_biting_rate * age_based_surface_area(age_in_years)

    # a = 1.933115996606731
    # b = 0.8970967578560122
    # c = 0.8991575730911449
    # d = -1.441447182704689

    # return 1 / (1 + np.exp(-(a * np.log(age_in_years) +
    #                          b * np.log(relative_biting_rate) +
    #                          c * np.log(daily_sim_eir) +
    #                          d)))
    a = 1.516316521864015
    b = 0.8930658890703256
    c = -0.042777966215856424

    return 1 / (1 + np.exp(-(a * np.log(age_in_years) +
                             b * np.log(individual_daily_eir) +
                             c)))

def infection_duration():
    pass
    # Calculate how long the person's infections last



if __name__ == "__main__":
    IntrahostComponent.set_params() # default params
    IntrahostComponent.update_params({"Run_Number": 2})


    ic = IntrahostComponent.create()
    # ic._configure_from_params({"Run_Number": 0})
    # ic.challenge()

    daily_sim_eir = 0.1


    print(f"1 challenge per {1/daily_sim_eir} days")

    all_data = []
    ab_conc = []
    ab_capacity = []
    pfemp_variant_frac = []
    num_infections = []
    total_pfemp_variants = 1070
    daily_eir_list = []
    infectiousness = []
    infection_stats = {}

    relative_biting_risk = 1

    final_age_years = 20

    for t in range(final_age_years*365):
        age_years = t/365
        surface_area = age_based_surface_area(age_years)

        # seasonal_effect = 0.45*np.sin(2*np.pi*t/365)+0.55
        seasonal_effect = 1
        daily_eir_with_seasonal_effect = daily_sim_eir*seasonal_effect

        # After 5 years, reduce the biting risk significantly
        # if age_years > 5:
        #     relative_biting_risk = 0.2

        individual_daily_eir = daily_eir_with_seasonal_effect * surface_area * relative_biting_risk

        bite_today = np.random.random() < individual_daily_eir

        if bite_today and ic.n_infections < 3:
            ic.challenge()

        # if t % int(1/individual_daily_eir) == 0:
        #     if ic.n_infections < 3:
        #         ic.challenge()

        # Keep track of each infection object and report when they are cleared
        infections_today = ic.infections
        infection_today_id_list = list(map(lambda x: id(x), infections_today))
        # print(infection_today_id_list)

        # If this is the first time an infection id has shown up in the list, add it to the infection_stats dictionary and record the time
        for infection_id in infection_today_id_list:
            if infection_id not in infection_stats:
                infection_stats[infection_id] = {"start_time": t}

        # If a previously seen infection id is not in the list today, record the end time, if it hasn't been recorded already
        for infection_id in infection_stats:
            if infection_id not in infection_today_id_list:
                if "end_time" not in infection_stats[infection_id]:
                    infection_stats[infection_id]["end_time"] = t

        ic.update(dt=1)

        n_active_this_timestep = len(ic.susceptibility.get_active_PfEMP1_major_antibodies())

        if n_active_this_timestep > 0:
            antibody_list = ic.susceptibility.get_active_PfEMP1_major_antibodies()
            antigen_count = list(map(lambda x: x.antigen_count, antibody_list))
            antibody_capacity = list(map(lambda x: x.antibody_capacity, antibody_list))
            antibody_concentration = list(map(lambda x: x.antibody_concentration, antibody_list))

            ab_conc.append(np.mean(antibody_concentration)*n_active_this_timestep/total_pfemp_variants)
            ab_capacity.append(np.mean(antibody_capacity)*n_active_this_timestep/total_pfemp_variants)
            pfemp_variant_frac.append(n_active_this_timestep/total_pfemp_variants)
            num_infections.append(ic.n_infections)

        else:
            ab_conc.append(np.nan)
            ab_capacity.append(np.nan)
            pfemp_variant_frac.append(0)
            num_infections.append(ic.n_infections)
        daily_eir_list.append(individual_daily_eir)
        infectiousness.append(ic.infectiousness)

        today_data = {"day": t,
                        "ab_conc": ab_conc[-1],
                        "ab_capacity": ab_capacity[-1],
                        "pfemp_variant_frac": pfemp_variant_frac[-1],
                        "num_infections": num_infections[-1],
                        "daily_eir": individual_daily_eir,
                        "infectiousness": ic.infectiousness}
        all_data.append(today_data)

    all_data_df = pd.DataFrame(all_data)

    pfemp1_variant_frac_predicted = predict_emod_pfemp1_variant_fraction(age_in_years=final_age_years,
                                                                         relative_biting_rate=1,
                                                                         daily_sim_eir=daily_sim_eir)
    print(f"Age in years = {age_years}")
    print("Predicted PfEMP1 variant fraction at t=3000: ", pfemp1_variant_frac_predicted)

    plt.figure()
    plt.subplot(321)
    plt.plot(ab_conc, label="Antibody Concentration")
    plt.plot(ab_capacity, label="Antibody Capacity")
    plt.xlabel("Time")
    plt.legend()

    plt.subplot(322)
    plt.plot(pfemp_variant_frac)
    plt.axhline(y=pfemp1_variant_frac_predicted, color='r', linestyle='--', label="Predicted PfEMP1 Variant Fraction")
    plt.xlabel("Time")
    plt.ylabel("Pfemp1VariantFraction")
    plt.legend()

    plt.subplot(323)
    plt.plot(num_infections, label="Number of infections")
    plt.xlabel("Time")
    plt.ylabel("num_infections")

    plt.subplot(324)
    plt.plot(daily_eir_list)
    plt.xlabel("Time")
    plt.ylabel("Daily EIR")

    plt.subplot(325)
    plt.plot(infectiousness)
    plt.xlabel("Time")
    plt.ylabel("Infectiousness")

    plt.tight_layout()
    plt.show()

    # print(infection_stats)

    plt.figure()
    for infection_id in infection_stats:
        # print(f"Infection {infection_id} started at {infection_stats[infection_id]['start_time']} and ended at {infection_stats[infection_id]['end_time']}")
        start_time = infection_stats[infection_id]['start_time']
        if 'end_time' in infection_stats[infection_id]:
            end_time = infection_stats[infection_id]['end_time']
            infection_stats[infection_id]['duration'] = end_time - start_time
            print(f"Duration: {end_time-start_time}")

            # Average infectiousness during this time
            avg_infectiousness = np.mean(infectiousness[start_time:end_time])
            print(f"Average infectiousness: {avg_infectiousness}")
            print(f"Peak infectiousness: {np.max(infectiousness[start_time:end_time])}")
            infection_stats[infection_id]['avg_infectiousness'] = avg_infectiousness
            infection_stats[infection_id]['peak_infectiousness'] = np.max(infectiousness[start_time:end_time])

            # How much of the total infectiousness comes from the peak?
            peak_infectiousness = np.max(infectiousness[start_time:end_time])
            total_infectiousness = np.sum(infectiousness[start_time:end_time])
            print(f"Peak infectiousness as a fraction of total infectiousness: {peak_infectiousness/total_infectiousness}")
            infection_stats[infection_id]['total_infectiousness'] = total_infectiousness


            plt.plot(infectiousness[start_time:end_time], c='black', alpha=0.1)
            plt.yscale("log")

    plt.show()


    infection_stats_df = pd.DataFrame.from_dict(infection_stats, orient='index')
    infection_stats_df.reset_index(inplace=True)
    infection_stats_df.columns = ['infection_id', 'start_time', 'end_time', 'duration', 'avg_infectiousness', 'peak_infectiousness', 'total_infectiousness']

    # Scatter plot of duration vs infectiousness, with points colored by order of appearance
    plt.figure()
    plt.scatter(infection_stats_df['duration'], infection_stats_df['avg_infectiousness'], c=range(len(infection_stats_df)))
    plt.xlabel('Duration')
    plt.ylabel('Average Infectiousness')
    plt.show()

