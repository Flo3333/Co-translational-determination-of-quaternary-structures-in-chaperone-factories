from bigfish.stack import check_parameter
import pandas as pd
import numpy as np
from CustomPandasFramework.operations import list_match, check_expectedcolumns, keep_columns

def clean_Sets_fromxl(file_fullpath) :
   """Prepare an Excel file for data explotation. Result is used for making TargetGenes and Oligopool Table.
   If Excel file format is incorrect raises an exception.
   
   Parameters
   ----------
      file_fullpath : str
         full path to excel file
         
   
   Returns
   -------
      res_Frame : pd.DataFrame
         Cleaned Frame ready for data exploitation
         
   """
   check_parameter(file_fullpath = (str))
   if file_fullpath[len(file_fullpath) - 4 : len(file_fullpath)] != "xlsx" :
      raise Exception("Entry file must be an Excel file with extension xlsx.")

   raw_Frame =  pd.read_excel(file_fullpath)
   expected_cols = ["ENSG", "ENST", "GeneName", "Set"]
   check_expectedcolumns(raw_Frame, expected_cols)

   #Cleaning
   res_Frame = keep_columns(raw_Frame, expected_cols)
   res_Frame = res_Frame.dropna(subset=["Set", "GeneName"])

   return res_Frame




def clean_TargetGenes(file_fullpath) :
   """Prepare an Excel file for data explotation. Result is used for making TargetGenes and Oligopool Table.
   If Excel file format is incorrect raises an exception.
   
   Parameters
   ----------
      file_fullpath : str
         full path to excel file
         
   
   Returns
   -------
      res_Frame : pd.DataFrame
         Cleaned Frame ready for data exploitation
         
   """

   check_parameter(file_fullpath = (str))
   if file_fullpath[len(file_fullpath) - 4 : len(file_fullpath)] != "xlsx" :
      raise Exception("Entry file must be an Excel file with extension xlsx.")
   
   raw_Frame =  pd.read_excel(file_fullpath, usecols= "A:D")
   res_Frame = raw_Frame.copy()

   if not list_match(list(raw_Frame.keys()), ["ENSG",	"ENST",	"GeneName",	"n° oligopool"]) :
      raise Exception("File placed in clean_TargetGenes has not expected format.")


   return res_Frame

def clean_PlatesPlan(file_fullpath) :
   """Prepare an Excel file for data explotation. Result is used for making Plates and list_PlatePlan.
   If Excel file format is incorrect raises an exception.
   
   Parameters
   ----------
      file_fullpath : str
         full path to excel file
         
   
   Returns
   -------
      res_Frame : pd.DataFrame
         Cleaned Frame ready for data exploitation
         
   """

   check_parameter(file_fullpath = (str))
   if file_fullpath[len(file_fullpath) - 4 : len(file_fullpath)] != "xlsx" :
      raise Exception("Entry file must be an Excel file with extension xlsx.")
   
   raw_Frame = pd.read_excel(file_fullpath, usecols= "D:R")
   raw_Frame.columns = ["D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"]
   clean_Frame = raw_Frame.drop("E", axis= 1)
   clean_Frame = clean_Frame.dropna(subset= "F").reset_index(drop= True)
   
   return clean_Frame

def list_PlatePlan(cleaned_PlatesPlan) :
   """Prepare a cleaned_PlatesPlan for Plates Dataframe making.

   Parameters
   ----------
      cleaned_PlatesPlan : pd.DataFrame
         DataFrame containing plates components informations : result from clean_PlatesPlan.
         
   
   Returns
   -------
      res_Frame : pd.DataFrame
         Frame containing the position of each plate plan
         
   """
   
   #Preparing dataframe
   if list(cleaned_PlatesPlan.keys()) != ["D","F","G","H","I","J","K","L","M","N","O","P","Q","R"] : 
      raise Exception("Unexpected data shape. Did you use clean_PlatesPlan?")

   PlatePlan_List = pd.DataFrame({
      "Oligopool Name" : [],
      "name" : [],
      "coord" : [] #expected tuple(x : axis0, y: axis1)
   })

   plates_start_idx = cleaned_PlatesPlan[cleaned_PlatesPlan["F"] == "A"].index
   Oligopool_name = []
   name = []
   coord = []

   for idx in plates_start_idx :
      if idx == 0 :
         new_name = np.NaN
      elif cleaned_PlatesPlan.at[idx-1, "F"] == "H" :
         new_name = np.NaN
      else :
         new_name = cleaned_PlatesPlan.at[idx-1, "F"]
        
      new_coord = (idx, "F")

      Oligopool_name += [cleaned_PlatesPlan.at[idx,"D"]]
      name += [new_name]
      coord += [new_coord]

   PlatePlan_List["Oligopool Name"] = Oligopool_name
   PlatePlan_List["name"] = name
   PlatePlan_List["coord"] = coord

   return PlatePlan_List




def convert_PlatePlan(cleaned_PlatesPlan, PlatePlan_coord):
   """Convert a PlatePlan to GeneMap which is needed for TargetGenes Frame making.

   Parameters
   ----------
      cleaned_PlatesPlan : pd.DataFrame
         DataFrame containing plates components informations : result from clean_PlatesPlan.
      PlatePlan_coord : tuple(int,str)
         coordinate of the plate plan withing the excel file.
         
   
   Returns
   -------
      res_Frame : pd.DataFrame
         Frame containing names and coordinates of each genes in Plate given by PlatePlan_coord.
         
   """
   


   mapping_col = {"G" : 1, "H" : 2, "I" : 3, "J": 4, "K":5, "L":6, "M":7 ,"N":8 ,"O":9, "P":10 ,"Q": 11, "R": 12}
   mapping_line = {1 : "A", 2 : "B", 3 : "C", 4 : "D", 5 : "E", 6 : "F", 7 : "G", 8 : "H"}

   GeneMap = pd.DataFrame({
   "GeneName" : [],
   "Coord" : [] # as expected in final data shape that is to say real coordinate during experimentation
   })

   Plate_Plan = cleaned_PlatesPlan.loc[PlatePlan_coord[0] : PlatePlan_coord[0] + 7, "G" : "R"].reset_index(drop= True)
   GeneName = []
   Coord = []

   for line in range(0,7) :
      for col in ["G","H","I","J","K","L","M","N","O","P","Q","R"]:
         GeneName += [Plate_Plan.at[line, col]]
         Coord += [(mapping_line[line + 1], mapping_col[col])]

   GeneMap["GeneName"] = GeneName
   GeneMap["Coord"] = Coord

   return GeneMap
