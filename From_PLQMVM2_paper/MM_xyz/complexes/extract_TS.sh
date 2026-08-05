#!/bin/bash
 lig_list='
 ac_67 ad_12 ad_17 ad_23 ad_24 ad_32
 ad_33 ad_63 ad_67 ad_70 ad_71 ad_73
 ad_74 ad_76 ad_78 ad_80 ad_81 ad_83
 kb_01 kb_02 kb_03 kb_04 kb_18 kb_19
 kb_53 kb_54 kb_55 kb_56 kb_57 kb_58
 kb_59 kb_60 kb_61 kb_62 kb_67 kb_69
 kk_98 kk_99'

 for lig in $lig_list; do
 ligname=umass_1_${lig}_2d_vm2_conformers.xyz
 name=2i0d_1580_p4a_tleap_${ligname}
 ts_info=complex_${lig}_TS_info.txt
 grep 'Rank=' $name > $ts_info
 done
