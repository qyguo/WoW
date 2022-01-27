As an example:

root -l -b -q slimNtuple.C\(2017,\"GluGluHToZZTo4L_M125_2017\",true,true,true,true,true\)

slimNtuple(const int & _year_=2017, const string & _name_DS_="GluGluHToZZTo4L_M125_2017", const bool & isMC = true,  const bool & isSignal = true, const bool & _OldNtuple=false, const bool & _addNewVar=false, const bool & _Test = false)

arg[0]: default year is 2017

arg[1]: default sample is "GluGluHToZZTo4L_M125_2017"

arg[2]: default isMC is true

arg[3]: default isSignal is true 

arg[4]: default _OldNtuple is false (_OldNtuple true will have less branch)

arg[5]: default _addNewVar is false (based on the existing old Ntuple to add new variables)

arg[6]: default _Test is false (just test 20k events)

