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

if __name__ == "__main__":
    IntrahostComponent.set_params() # default params
    IntrahostComponent.update_params({"Run_Number": 4})

    N_humans_to_simulate = 100
    daily_sim_eir = 0.10
    final_age_years = 40
    print(f"1 challenge per {1/daily_sim_eir} days")

    human_infection_df_list = []

    for j in range(N_humans_to_simulate):
        print("Human id = ", j)
        ic = IntrahostComponent.create()

        all_data = []
        pfemp_variant_frac = []
        num_infections = []
        total_pfemp_variants = 1070
        daily_eir_list = []
        infection_stats = {}

        relative_biting_risk = 1
        gametocyte_densities = {}

        custom_infection_id = 0
        infection_ids_to_pop = []
        for t in range(final_age_years*365):
            age_years = t/365
            if age_years % 5 == 0:
                print(f"Age in years: {age_years}")
            surface_area = age_based_surface_area(age_years)

            daily_eir_with_seasonal_effect = daily_sim_eir
            individual_daily_eir = daily_eir_with_seasonal_effect * surface_area * relative_biting_risk

            bite_today = np.random.random() < individual_daily_eir

            if bite_today and ic.n_infections < 3:
                ic.challenge()
                # print("challenged")

            # Keep track of each infection object and report when they are cleared
            infections_today = ic.infections
            infection_today_id_list = list(map(lambda x: id(x), infections_today))
            gametocyte_densities_today = list(map(lambda x: x.mature_gametocyte_density, infections_today))

            # If this is the first time an infection id has shown up in the list, add it to the infection_stats dictionary and record the time
            for infection_id in infection_today_id_list:
                if infection_id not in infection_stats:
                    infection_stats[infection_id] = {"start_time": t}
                    gametocyte_densities[infection_id] = []

            # Keep track of gametocyte density for each current infection:
            for iid, g in zip(infection_today_id_list, gametocyte_densities_today):
                gametocyte_densities[iid].append(g)

            # If a previously seen infection id is not in the list today, record the end time, if it hasn't been recorded already
            for infection_id in infection_stats:
                if infection_id not in infection_today_id_list:
                    if "end_time" not in infection_stats[infection_id]:
                        infection_stats[infection_id]["end_time"] = t
                        infection_stats[infection_id]["duration"] = t - infection_stats[infection_id]["start_time"]
                        gametocyte_densities[infection_id] = sum(gametocyte_densities[infection_id])

                        infection_ids_to_pop.append(infection_id)

            # Rename infection ID to avoid collision
            for infection_id in infection_ids_to_pop:
                infection_stats[custom_infection_id] = infection_stats.pop(infection_id)
                gametocyte_densities[custom_infection_id] = gametocyte_densities.pop(infection_id)
                custom_infection_id += 1
            infection_ids_to_pop = []



            ic.update(dt=1)

            n_active_this_timestep = len(ic.susceptibility.get_active_PfEMP1_major_antibodies())

            if n_active_this_timestep > 0:
                antibody_list = ic.susceptibility.get_active_PfEMP1_major_antibodies()
                antigen_count = list(map(lambda x: x.antigen_count, antibody_list))
                antibody_capacity = list(map(lambda x: x.antibody_capacity, antibody_list))
                antibody_concentration = list(map(lambda x: x.antibody_concentration, antibody_list))

                pfemp_variant_frac.append(n_active_this_timestep/total_pfemp_variants)
                num_infections.append(ic.n_infections)

            else:
                pfemp_variant_frac.append(0)
                num_infections.append(ic.n_infections)
            daily_eir_list.append(individual_daily_eir)

            today_data = {"day": t,
                          "pfemp_variant_frac": pfemp_variant_frac[-1],
                          "num_infections": num_infections[-1],
                          "daily_eir": individual_daily_eir,
                          "infectiousness": ic.infectiousness}
            all_data.append(today_data)

        all_data_df = pd.DataFrame(all_data)

        # Ignore any current infections
        for infection_id in infection_today_id_list:
            del infection_stats[infection_id]
            del gametocyte_densities[infection_id]

        # For each infection, get duration and total gametocyte density:
        iids = []
        durations = []
        total_gametocyte_densities = []
        start_days = []
        for infection_id in infection_stats:
            iids.append(infection_id)
            start_days.append(infection_stats[infection_id]['start_time'])
            durations.append(infection_stats[infection_id]['end_time'] - infection_stats[infection_id]['start_time'])
            total_gametocyte_densities.append(gametocyte_densities[infection_id])

        df_infections = pd.DataFrame({"infection_id": iids,
                                      "start_day": start_days,
                                      "duration": durations,
                                      "total_gametocyte_density": total_gametocyte_densities})

        # Merge all_data_df and df_infections to get the pfemp_variant_frac at the start time of each infection
        df_infections = pd.merge(all_data_df[["day", "pfemp_variant_frac"]], df_infections, left_on="day", right_on="start_day", how="right")

        df_infections["human_id"] = j
        human_infection_df_list.append(df_infections)

    df_all = pd.concat(human_infection_df_list)

    print("Saving infection data to infection_data.csv")
    df_all.to_csv("infection_data.csv", index=False)

    plt.figure()
    plt.subplot(221)
    plt.scatter(df_all['start_day'], df_all['pfemp_variant_frac'], c=df_all["human_id"], alpha=0.1)
    plt.xlabel("Infection start day")
    plt.ylabel("PfEMP1 Variant Fraction")
    plt.subplot(222)
    plt.scatter(df_all['pfemp_variant_frac'], df_all['duration'], c=df_all["human_id"], alpha=0.1)
    plt.xlabel("PfEMP1 Variant Fraction")
    plt.ylabel("Infection Duration")
    plt.subplot(223)
    plt.scatter(df_all['pfemp_variant_frac'], df_all['total_gametocyte_density'], c=df_all["human_id"], alpha=0.1)
    plt.xlabel("PfEMP1 Variant Fraction")
    plt.ylabel("Total Gametocyte Density")
    plt.subplot(224)
    plt.scatter(df_all['pfemp_variant_frac'], df_all['total_gametocyte_density']/df_all['duration'], c=df_all["human_id"], alpha=0.1)
    plt.xlabel("PfEMP1 Variant Fraction")
    plt.ylabel("Mean Gametocyte Density per Day")
    plt.tight_layout()
    plt.show()
