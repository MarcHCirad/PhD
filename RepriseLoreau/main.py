import numpy as np
import matplotlib.pyplot as plt


def recruitment(Ressource, beta, a, c):
    rslt = beta * np.exp(-(Ressource - a)**2 / (c*Ressource)) / Ressource**1.5
    return rslt

def death(Ressource, delta, delta_min, rn):
    return delta / (1+ np.exp(Ressource-rn)) + delta_min

def main():
    # rn, delta, delta_min = 4, 0.0535, 0.004
    beta, a, c = 0.04, 10, 6

    # Ressource_tab = np.arange(0.1, 12, 0.1)
    # Death = np.array([death(R, delta, delta_min, rn) for R in Ressource_tab])
    # Recruitment = np.array([recruitment(R, beta, a, c) for R in Ressource_tab])
    # Growth_rate = Recruitment-Death
    # print(max(Recruitment), max(Growth_rate))
    # plt.plot(Ressource_tab, Growth_rate, Ressource_tab, Recruitment, Ressource_tab, Death)
    # plt.legend(["Growth", "Recruitment", "Death"])
    # plt.show()
    print(recruitment(5, beta, a, c))


if __name__ == "__main__":
    main()