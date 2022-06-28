nohup python -u extract_JES_nuis_2D.py --obsName="njets_pt30_eta4p7 vs pT4l" --obsBins="|0|1| vs |0|15|30|13000| / |1|2| vs |0|60|80|120|13000| / |2|10| vs |0|100|170|250|13000|"  --year="2016" >& log_extract_JES_njets_pt30_eta4p7_vs_pT4l_2016.txt &  # 5 bins

#nohup python -u extract_JES_nuis_2D.py --obsName="pt_leadingjet_pt30_eta4p7 vs pTj2" --obsBins="|-2.0|30.0| vs |-11|30| / |30|60| vs |30|60| / |60|350| vs |30|60| / |60|350| vs |60|350|"  --year="2016" >& log_extract_JES_pt_leadingjet_pt30_eta4p7_vs_pTj2_2016.txt &  # 5 bins


nohup python -u extract_JES_nuis_2D.py --obsName="TauC_Inc_0j_EnergyWgt vs pT4l" --obsBins="|-1|15| vs |0|15|30|45|70|120|13000| / |15|25| vs |0|120|13000| / |25|40| vs |0|120|13000| / |40|13000| vs |0|200|13000|"  --year="2016" >& log_extract_JES_TauC_Inc_0j_EnergyWgt_vs_pT4l_2016.txt &  # 5 bins


#nohup python -u extract_JES_nuis_2D.py --obsName="pT4l vs pT4lj" --obsBins="|0|13000| vs |-2|0| / |0|85| vs |0|30| / |85|350| vs |0|45| / |0|85| vs |30|350| / |85|350| vs |45|350|"  --year="2016" >& log_extract_pT4l_vs_pT4lj_2016.txt &  # 5 bins


