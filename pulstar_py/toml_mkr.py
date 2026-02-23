### This will be my working code to create toml strings that will be parsed into the 
##rust program. 
from pydantic import BaseModel,Field
from typing import Union,Optional,List
import polars as pl

#----------------------------------------
#----------Pulstar input-----------------
#----------------------------------------

#define pulsation mode
class Mode(BaseModel):
    l:int
    m:int
    rel_dr:float
    k:float
    frequency:float
    phase_offset:float
    rel_dtemp:float
    phase_rel_dtemp:float
    rel_dg:float
    phase_rel_dg:float
#    def __init__(self,l=4,m=1,rel_dr=0.024,k=0.05,frequency=6.74,phase_offset=0.0,rel_dtemp=2.62,phase_rel_dtemp=180.0,rel_dg=10.0,phase_rel_dg=34.0):
#        self.l=l
#        self.m=m
#        self.rel_dr=rel_dr
#        self.k=k
#        self.frequency=frequency
#        self.phase_offset=phase_offset
#        self.rel_dtemp=rel_dtemp
#        self.phase_rel_dtemp=phase_rel_dtemp
#        self.rel_dg=rel_dg
#        self.phase_rel_dg=phase_rel_dg

#define star data
class StarData(BaseModel):
    mass:float
    radius:float
    effective_temperature:float
    v_omega:float
    inclination_angle:float
#    def __init__(self, mass=10.0,radius=6.93,effective_temperature=22000.0,v_omega=20.0,inclination_angle=45.0):
#        self.mass=mass
#        self.radius=radius
#        self.effective_temperature=effective_temperature
#        self.v_omega=v_omega
#        self.inclination_angle=inclination_angle
    
#define phases of pulsation
class UniformTime(BaseModel):
    start:float
    end:float
    step:float
#    def __init__(self, start=0.0,end=0.1,step=0.01):
#        self.start=start
#        self.end=end
#        self.step=step

class ExplicitTime(BaseModel):
    collection:List[float]
#    def __init__(self,collection=[0.01,0.3,0.32]):
#        self.collection=collection

class TimePoints(BaseModel):
    Uniform:Optional[UniformTime]=None
    Explicit:Optional[ExplicitTime]=None
#    def __init__(self,Uniform=UniformTime(start=0.0,end=0.1,step=0.01),Explicit=None):
#        self.Uniform=Uniform
#        self.Explicit=Explicit

#define type of meshing
class SphericalStar(BaseModel):
    theta_step:float
    phi_step:float
#    def __init__(self,theta_step=2.0,phi_step=4.0):
#        self.theta_step=theta_step
#        self.phi_step=phi_step

class Mesh(BaseModel):
    Sphere:Optional[SphericalStar]=None
#    def __init__(self,sphere=SphericalStar(theta_step=2.0,phi_step=4.0)):
#        self.Sphere=sphere



class PulstarConfig(BaseModel):
    mode_data:List[Mode]
    star_data:StarData
    time_points:TimePoints
    mesh:Mesh
#    def __init__(self,mode_data=[(Mode())],star_data=StarData(),time_points=TimePoints(),mesh=Mesh()):
#        self.mode_data=mode_data
#        self.star_data=star_data
#        self.time_points=time_points
#        self.mesh=mesh

#----------------------------------------
#----------Profile input-----------------
#----------------------------------------

#Intensity grid
class JorisGrid(BaseModel):
    temperature:float
    log_gravity:float
    filename:str

class NadyaGrid(BaseModel):
    temperature:float
    log_gravity:float
    metalicity:float
    filename:str


class IntensityGrid(BaseModel):
    Joris:Optional[JorisGrid]=None
    Nadya:Optional[NadyaGrid]=None

#Wave length range
class WavelengthRange(BaseModel):
    start:float
    end:float
    step:float

#profile_config
class ProfileConfig(BaseModel):
    max_velocity:float
    path_to_grids:str
    wavelength_range:WavelengthRange
    intensity_grids:List[IntensityGrid]


import tomli_w

#mode=Mode(l=4,m=1,rel_dr=0.024,k=0.05,frequency=6.94,phase_offset=0.03,rel_dtemp=)
mode=Mode(l=4,m=1,rel_dr=0.024,k=0.05,frequency=6.74,phase_offset=0.0,rel_dtemp=2.62,phase_rel_dtemp=180.0,rel_dg=10.0,phase_rel_dg=34.0)
star_data=StarData(mass=10.0, radius=6.94, effective_temperature=22000.0,v_omega=20.0,inclination_angle=45.0)
time_points=TimePoints(Uniform=UniformTime(start=0.0,end=0.1,step=0.01))
mesh=Mesh(Sphere=SphericalStar(theta_step=2.0,phi_step=4.0))

star_data.radius=6.93
pulsconfig=PulstarConfig(mode_data=[mode],star_data=star_data,time_points=time_points,mesh=mesh)

pulsconfig_dict=pulsconfig.model_dump(exclude_none=True)

puls_toml_string=tomli_w.dumps(pulsconfig_dict)


####Profile config
wl_range=WavelengthRange(start=413.85,end=414.0,step=0.002)
path_to_grids="../profile/grids/"
grid1=IntensityGrid(Joris=JorisGrid(temperature=21000.0,log_gravity=3.5,filename="t21000g35.txt"),Nadya=None)
grid2=IntensityGrid(Joris=JorisGrid(temperature=21000.0,log_gravity=4.5,filename="t21000g45.txt"),Nadya=None)
grid3=IntensityGrid(Joris=JorisGrid(temperature=24000.0,log_gravity=3.5,filename="t24000g35.txt"),Nadya=None)
grid4=IntensityGrid(Joris=JorisGrid(temperature=24000.0,log_gravity=4.5,filename="t24000g45.txt"),Nadya=None)
prof_config=ProfileConfig(max_velocity=1.0e2,path_to_grids=path_to_grids,wavelength_range=wl_range,intensity_grids=[grid1,grid2,grid3,grid4])

prof_config_dict=prof_config.model_dump(exclude_none=True)

prof_toml_string=tomli_w.dumps(prof_config_dict)

#print(puls_toml_string)
#print("and then")

#print(prof_toml_string)


print("-----------")
print("-----------")
print("now for the real test")

print("-----------")
print("-----------")

import pulstar_py as pls_py
pls_py.propulse(prof_toml_string,puls_toml_string)

df = pl.read_parquet("wavelengths_tp10.parquet")
df_old = pl.read_parquet("wave_old.parquet")
print(df.head(5))
print(df.head(5))
