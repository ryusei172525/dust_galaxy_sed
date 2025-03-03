# Dust to gas mass ratio 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import astropy.units as u
import pandas as pd
import matplotlib.pyplot as plt
import os
from IPython.display import display

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

# TODO: 観測データ表示はデザインの見直しが必要
def observational_sed_diagram(observational_data):
    # データの読み込み
    data = np.loadtxt(observational_data, skiprows=1)
    t = data[:, 0]
    y = data[:, 1]
    error_plus = data[:, 2]
    error_minus = data[:, 3]
    upper_limit = data[:, 4].astype(bool)

    # 通常のデータと上限データを分ける
    t_normal = t[~upper_limit]
    y_normal = y[~upper_limit]
    errors_normal = [error_minus[~upper_limit], error_plus[~upper_limit]]

    t_upper = t[upper_limit]
    y_upper = y[upper_limit]

    # プロットの設定
    plt.figure(figsize=(10, 6))

    # 通常の観測データをプロット（エラーバー付き）
    plt.errorbar(t_normal, y_normal, yerr=errors_normal, fmt='o', color='blue', label='Observed Data')

    # 上限データに矢印と水平線を追加
    plt.scatter(t_upper, y_upper, marker='v', color='red', label='Upper Limit')
    for i in range(len(t_upper)):
        # 水平線の追加
        plt.hlines(y_upper[i], t_upper[i] - 0.1 * t_upper[i], t_upper[i] + 0.1 * t_upper[i], color='red', linestyles='dashed')
        
        # 矢印の追加（矢印先端が水平線の下に来るよう調整）
        plt.arrow(t_upper[i], y_upper[i] * 0.8, 0, -0.3 * y_upper[i], head_width=0.02 * t_upper[i], 
                  head_length=0.1 * y_upper[i], fc='red', ec='red')

    # グラフの設定
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Wavelength [$\mu$m]', fontsize=15)
    plt.ylabel('$\lambda L_\lambda $[erg s$^{-1}$ $L_\odot^{-1}$]', fontsize=15)
    plt.title('Observation Data with Upper Limits and Error Bars')
    plt.legend()
    plt.grid()
    plt.xlim(0.2e-1, 1e3)
    plt.ylim(1e2, 1e13)

    plt.show()

def sed_diagram_with_observational_data(observational_data):
    # 図のサイズを横長に設定
    plt.figure(figsize=(12, 6))

    # 既存の SED 図をプロット
    sed_diagram()

    # 観測データをロード
    data = np.loadtxt(observational_data, skiprows=1)
    t = data[:, 0]  # 波長
    y = data[:, 1]  # 明るさ
    error_plus = data[:, 2]
    error_minus = data[:, 3]
    upper_limit = data[:, 4].astype(bool)

    # 通常データと上限データに分ける
    t_normal = t[~upper_limit]
    y_normal = y[~upper_limit]
    errors_normal = [error_minus[~upper_limit], error_plus[~upper_limit]]

    t_upper = t[upper_limit]
    y_upper = y[upper_limit]

    # 既存のプロットに追加
    ax = plt.gca()  # 現在の図の軸を取得

    # 通常の観測データをエラーバー付きでプロット
    ax.errorbar(t_normal, y_normal, yerr=errors_normal, fmt='o', color='blue', label='Observed Data')

    # 上限データをプロット（破線と下向き矢印）
    for i in range(len(t_upper)):
        ax.hlines(y_upper[i], t_upper[i] * 0.9, t_upper[i] * 1.1, color='red', linestyles='dashed')
        ax.annotate('', xy=(t_upper[i], y_upper[i] * 0.5), xytext=(t_upper[i], y_upper[i] * 12),
                    arrowprops=dict(arrowstyle="->", color='red'))

    # 軸ラベルと設定の追加
    ax.set_xlabel('Wavelength [$\mu$m]', fontsize=15)
    ax.set_ylabel('$\lambda L_\lambda $[erg s$^{-1}$ $L_\odot^{-1}$]', fontsize=15)

    # 凡例を図の外に配置
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), borderaxespad=0., fontsize=12)

    # 軸スケールと範囲
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.2e-1, 1e3)
    ax.set_ylim(1e6, 1e13)

    plt.tight_layout()  # レイアウト調整
    plt.show()

def find_file(directory, pattern):
    files = []
    for root, dirs, filenames in os.walk(directory):
        for filename in filenames:
            if filename.startswith(pattern):
                files.append(os.path.join(root, filename))
    return files
    
def MassDataTable(galaxy_mass):
    directory = "galaxy_sed/Dust_cpp_pybind/calculation_result/dust_model/asano_model/"
    pattern = "total_dust_mass_t_sf"
    file = find_file(directory, pattern)
    
    # Check if exactly one file is found
    if len(file) == 1:
        file_path = file[0]
        
        df = pd.read_csv(file_path, delim_whitespace=True, header=0)

        galaxy_mass = galaxy_mass * 1e9
        # Modify columns based on galaxy_mass
        # df["M_galaxy[Msun]"] = galaxy_mass
        # df["M_gas[Msun]"] *= galaxy_mass
        # df["M_star[Msun]"] *= galaxy_mass
        # df["m_dust[Msun]"] *= galaxy_mass
        
        display(df)
        
        # 指定したカラムの値を行ごとに足し合わせる
        # df['sum_values'] = df['M_gas[Msun]'] + df['M_star[Msun]'] + df['m_gas_carbon[Msun]'] + df['m_gas_silicate[Msun]'] + df['M_Z_val[Msun]'] + df['m_dust_C[Msun]'] + df['m_dust_Si[Msun],'] + df['m_dust_Fe[Msun],'] + df['m_dust[Msun]']
        df['M_gas + M_star'] = df['M_gas[Msun]'] + df['M_star[Msun]']
        df['gas: carbon + silicate'] = df['m_gas_carbon[Msun]'] + df['m_gas_silicate[Msun]']
        df['dust: C + Si + Fe'] = df['m_dust_C[Msun]'] + df['m_dust_Si[Msun],'] + df['m_dust_Fe[Msun],']


        # 結果の確認
        print("df['M_gas + M_star']\n",df['M_gas + M_star'])
        print("df['gas: carbon + silicate']\n",df['gas: carbon + silicate'])
        print("df['dust: C + Si + Fe']\n",df['dust: C + Si + Fe'])
        
    else:
        print("Error: Expected exactly one matching file, but found", len(file))
