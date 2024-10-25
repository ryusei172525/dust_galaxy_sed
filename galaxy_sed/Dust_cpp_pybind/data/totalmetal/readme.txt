totalmetalA,B,Cディレクトリに入っているデータは

AGBstars van den Hoek & Groenewegen (1997)
mass range 1–7 M and metallicities Z = (5.0*10^-2,0.2,0.4,1.0) Zsun

SNe II Woosley & Weaver (1995)
mass range 12–40 M and metallicities Z = (5.0*10^-2,0.1,1.0) Zsun.
これらのデータを線型補間(内挿)してデータを作成している。

多分うえのふたつか関係ない（2020/11/25 NISHIDA)

AGBmetaldata.odsはKarakas 2010のデータである。
左から、
列の数、
星質量(Msun)、
Partial mixing zones: yes (in original table A6) or no、
metallicity、
final mass(Msun)、
species i、
Atomic mass of the species i、
Net yield(Msun)、
mass of species i lost in the wind(Msun)、
mass of species i initially present in the wind(Msun)、
average abundance of species in the wind (mass fraction) <X(i)>、
initial mass fraction of species i (mass fraction) X0(i)、
priduction factor:log10(<X(i)>/X0(i))
このうち、mass of species i lost in the windのデータを
metalの供給量であると、Inoue 2011では計算している。

AGBmetaldata_z**.datはKarakas(2010)の各metallicityのAGBデータである。
左から上記の列の数を除いたものである。

SNmetaldataz**.odsはNomoto et al.(2006)のデータである。
最上段は星質量(Msun)
次段は爆発エネルギー(E51)
下段はremnant mass(Msun)
以下は各speciesの放出量(Msun)
また、fileのz以下の数字はmetallicityを表している。

これらはInoue 2011で使用されたデータである。
各元素を見ていくと、AGBでは作らず、SNで作るという
元素が散見する。

datafile1.txtはKobayashi et al.(2006)のデータである。
各列の説明はファイルの中に記述してある。
M_finalはSNが起こる前の質量であり、progenitor massではないことに注意。
progenitor massと違う理由はstellar windによって吹き飛ばされたためである。
M_cutはremnant massである。
なので、放出されるmassはprogenitor mass - M_cutである。
放出されるmetal massはM_final - M_cut - M_H - M_Heである。
stellar windには金属が入っていないかどうかについてだが、
放出されるのは主にH, He layerなので金属量には依存しない？
それとも放出の部分で加わっている？

SNmetal_kobayashi_z**.datはKobayashi et al.(2006)の
各metallicityのデータである。
左からのデータはdatafile1.txtを参照

newyieldディレクトリに入っているデータは
Karakas (2010) : AGB_DATA_TYPE, Kobayashi et al.(2006) : SNeの
データを線型補間(内挿)してデータを作成している。
stellar metallicityの値などはAGBmetaldata.ods(AGB_DATA_TYPE)、datafile1.txt(SN)を
参照のこと。

Karakas (2010)とKobayashi et al.(2006)のデータ補間の手順
1. interpolate_metal_**.datを用いてそれぞれのデータをmetallicity方向に線型補間
　 補間したデータはnewyield/**metal_interpolate/metal_interpolateに入っている)
   なお、このディレクトリに入っているデータは左から星質量(Msun),mass of elements more heavier than He (Msun),
   Carbon mass (Msun), Silicon mass (Msun)となっている。

2. newyield/**metal_interpolate/metal_interpolateに入っているデータを用いて
   星質量方向に0.01Msun刻みで線型補間。ただし、6-13MsunはそれぞれAGB,SNのデータを用いて線型補間する。
   metallicityに関しては、AGBが0.005Zsunからしかないので、SNもこれに合わせる。
   使用したプログラムはinterpolate_smass.c
