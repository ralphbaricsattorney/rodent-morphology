import polars as pl
import numpy as np


def extract_std(aggregate_record):
    aggregate_record = aggregate_record.split("±")
    return float(aggregate_record[1].replace(" ", "").replace("*", "").replace("g", ""))


def extract_mean(aggregate_record):
    aggregate_record = aggregate_record.split("±")
    return float(aggregate_record[0].replace(" ", "").replace("*", ""))


def calculate_effect_sizes(
    male_n,
    female_n,
    male_means,
    female_means,
    male_variances,
    female_variances,
):
    cohen_ds = male_means - female_means
    pooled_sd = (
        ((male_n - 1) * male_variances) + ((female_n - 1) * female_variances)
    ) / (male_n + female_n - 2)
    try:
        cohen_ds /= pl.DataFrame(np.sqrt(pooled_sd))
    except:
        cohen_ds = pl.DataFrame(cohen_ds)
        cohen_ds / pl.DataFrame(np.sqrt(pooled_sd))
    return cohen_ds


def process_house_mouse_cranial_data():
    data = pl.read_csv("data/house_mouse_morphology.csv")
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[2])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[4])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[2])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[4])).rename({"map": "female_std"})
    )
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    male_n = data.select(["Males (N)"]).transpose(include_header=True)[:, 1:]
    female_n = data.select(["Females (N)"]).transpose(include_header=True)[:, 1:]
    effect_sizes = calculate_effect_sizes(
        male_n, female_n, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_guinea_pig_volume_data():
    data = pl.read_csv("data/guinea_pig_volume_morphology.csv")[:, :5]
    male_variances = pl.DataFrame(
        np.square(data.select(["male_volume_SD"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(
            data.select(["female_volume_SD"]).transpose(include_header=True)[:, 1:]
        )
    )
    male_means = data.select(["male_volume_Mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_volume_Mean"]).transpose(include_header=True)[
        :, 1:
    ]
    effect_sizes = calculate_effect_sizes(
        8, 8, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_guinea_pig_neuronal_data():
    data = pl.read_csv("data/guinea_pig_neuron_morphology.csv")[:, :5]
    male_variances = pl.DataFrame(
        np.square(
            data.select(["male_neuron_density_SD"]).transpose(include_header=True)[
                :, 1:
            ]
        )
    )
    female_variances = pl.DataFrame(
        np.square(
            data.select(["female_neuron_density_SD"]).transpose(include_header=True)[
                :, 1:
            ]
        )
    )
    male_means = data.select(["male_neuron_density_Mean"]).transpose(
        include_header=True
    )[:, 1:]
    female_means = data.select(["female_neuron_density_Mean"]).transpose(
        include_header=True
    )[:, 1:]
    effect_sizes = calculate_effect_sizes(
        8, 8, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_black_rat_data():
    data = pl.read_csv("data/rattus_rattus_morphology.csv")[:, :-1]
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[1])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[0])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[1])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[0])).rename({"map": "female_std"})
    )
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    effect_sizes = calculate_effect_sizes(
        14, 10, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_brown_rat_data():
    data = pl.read_csv("data/rattus_norvegicus_morphology.csv")[1:, :]
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[2])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[1])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[2])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[1])).rename({"map": "female_std"})
    )
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    effect_sizes = calculate_effect_sizes(
        58, 62, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_nutria_data():
    data = pl.read_csv("data/nutria_intestinal_morphology.csv", columns=[0, 1, 2])
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[1])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[2])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[1])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[2])).rename({"map": "female_std"})
    )
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    effect_sizes = calculate_effect_sizes(
        12, 18, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_mongolian_hamster_data():
    data = pl.read_csv("data/mongolian_hamster_morphology.csv")
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[1])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[2])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[1])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[2])).rename({"map": "female_std"})
    )
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    effect_sizes = calculate_effect_sizes(
        49, 55, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_mongolian_gerbil_data():
    data = pl.read_csv("data/mongolian_gerbil_morphology.csv")[1:, :3]
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[1])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[2])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[1])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[2])).rename({"map": "female_std"})
    )
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    effect_sizes = calculate_effect_sizes(
        61, 84, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_guinea_pig_data():
    data = pl.read_csv("data/guinea_pig_morphology.csv", columns=[0, 1, 2])[0, :]
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[1])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[2])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[1])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[2])).rename({"map": "female_std"})
    )
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    effect_sizes = calculate_effect_sizes(
        15, 17, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_golden_hamster_data():
    data = pl.read_csv("data/golden_hamster_morphology.csv", columns=[0, 3, 4])[1:, :]
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[1])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[2])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[1])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[2])).rename({"map": "female_std"})
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    effect_sizes = calculate_effect_sizes(
        65, 25, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_chinchilla_data():
    data = pl.read_csv("data/chinchilla_cranial_morphology.csv", columns=[0, 1, 3])
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[1])).rename({"map": "male_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_mean(s[2])).rename({"map": "female_mean"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[1])).rename({"map": "male_std"})
    )
    data = data.with_columns(
        data.map_rows(lambda s: extract_std(s[2])).rename({"map": "female_std"})
    )
    male_means = data.select(["male_mean"]).transpose(include_header=True)[:, 1:]
    female_means = data.select(["female_mean"]).transpose(include_header=True)[:, 1:]
    male_variances = pl.DataFrame(
        np.square(data.select(["male_std"]).transpose(include_header=True)[:, 1:])
    )
    female_variances = pl.DataFrame(
        np.square(data.select(["female_std"]).transpose(include_header=True)[:, 1:])
    )
    effect_sizes = calculate_effect_sizes(
        49, 47, male_means, female_means, male_variances, female_variances
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)
    return effect_sizes, variability_ratios


def process_black_tailed_prairie_dog_data():
    data = pl.read_csv("data/black_tailed_prairie_dog_morphology.csv")
    males = data.filter(pl.col("SEX") == "M")
    females = data.filter(pl.col("SEX") == "F")

    male_means = males.mean()[:, -1]
    female_means = females.mean()[:, -1]

    male_variances = males.var()[:, -1]
    female_variances = females.var()[:, -1]

    effect_sizes = calculate_effect_sizes(
        len(males),
        len(females),
        male_means,
        female_means,
        male_variances,
        female_variances,
    )
    variability_ratios = pl.DataFrame(male_variances / female_variances)

    return effect_sizes, variability_ratios


def process_ansells_mole_rat_data():
    data = pl.read_csv("data/ansells_mole_rat_cranial_morphology.csv")
    males = data.filter(pl.col("Sex") == "m")
    females = data.filter(pl.col("Sex") == "f")

    male_means = males.mean()[:, 4:]
    female_means = females.mean()[:, 4:]

    male_variances = males.var()[:, 4:]
    female_variances = females.var()[:, 4:]

    effect_sizes = calculate_effect_sizes(
        len(males),
        len(females),
        male_means,
        female_means,
        male_variances,
        female_variances,
    )
    variability_ratios = male_variances / female_variances

    return effect_sizes, variability_ratios
