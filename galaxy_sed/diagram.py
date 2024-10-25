# Dust to gas mass ratio 

import pandas as pd
import matplotlib.pyplot as plt
import glob
import astropy.units as u
import pandas as pd
import matplotlib.pyplot as plt
import os

sun_luminosity = 3.828e33 * u.erg/u.s

def delete_past_data():
    path = "galaxy_sed/Dust_cpp_pybind/calculation_result/"
    
    # 指定されたパス内のすべてのファイルとサブディレクトリを走査
    for root, dirs, files in os.walk(path):
        # ファイルを削除
        for file in files:
            file_path = os.path.join(root, file)
            try:
                os.remove(file_path)  # ファイルを削除
                print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Error deleting {file_path}: {e}")

def dust_to_gas_mass_ratio():
    path = "galaxy_sed/Dust_cpp_pybind/calculation_result/dust_model/asano_model/total*"
    file = glob.glob(path)[0]

    if not file:
        print("No file found matching the pattern.")
        return  # もしファイルが見つからなかったら、関数を終了

    dataframe = pd.read_csv(file, sep=' ')

    plt.plot(dataframe['age[Myr]'], dataframe['m_dust[Msun]'] / dataframe['M_gas[Msun]'], linestyle="solid", color="blue")

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Galaxy age (Myr)')
    plt.ylabel('Dust to gas mass ratio')
    plt.show()  # プロットを表示

def sed_diagram():
    path="galaxy_sed/Dust_cpp_pybind/calculation_result/makeSED/SED*"
    file = glob.glob(path)[0]

    if not file:
        print("No file found matching the pattern.")
        return  # もしファイルが見つからなかったら、関数を終了

    dataframe = pd.read_csv(file, sep=' ')

    plt.plot(dataframe['Wavelength(cm)']*1e4, dataframe['Total']*dataframe['Wavelength(cm)']/sun_luminosity, label='Total(stellar + dust grains SED)',linestyle = "solid", color="black")
    plt.plot(dataframe['Wavelength(cm)']*1e4, dataframe['Stellar_int']*dataframe['Wavelength(cm)']/sun_luminosity, label='intrinsic stellar SED',linestyle = "solid", color="blue")
    plt.plot(dataframe['Wavelength(cm)']*1e4, dataframe['Sil']*dataframe['Wavelength(cm)']/sun_luminosity, label='Silicate',linestyle = "solid", color="orange")
    plt.plot(dataframe['Wavelength(cm)']*1e4, dataframe['Gra']*dataframe['Wavelength(cm)']/sun_luminosity, label='Graphite',linestyle = "solid", color="green")
    plt.plot(dataframe['Wavelength(cm)']*1e4, dataframe['PAHneu']*dataframe['Wavelength(cm)']/sun_luminosity, label='neutral PAH',linestyle = "solid", color="purple")
    plt.plot(dataframe['Wavelength(cm)']*1e4, dataframe['PAHion']*dataframe['Wavelength(cm)']/sun_luminosity, label='ionized PAH',linestyle = "solid", color="red")

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Wavelength [$\mu$m]', fontsize=15)
    plt.ylabel('$\lambda L_\lambda $[erg s$^{-1}$ $L_\odot^{-1}$]', fontsize=15)
    # 凡例を図の中に配置
    plt.legend(loc='upper left', bbox_to_anchor=(0.35, 0.42))
    # plt.legend(bbox_to_anchor=(1.05, 1.0),loc='upper left', fontsize=9)
    plt.xlim(0.8e-1,1e3)
    plt.ylim(1e1,1e13)