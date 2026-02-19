### This will be my working code to create toml strings that will be parsed into the 
##rust program. 
from pydantic import BaseModel, RootModel
from typing import Union,Optional,List

#define pulsation mode
class mode(BaseModel):
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
    def __init__(self,l=4,m=1,rel_dr=0.024,k=0.05,frequency=6.74,phase_offset=0.0,rel_dtemp=2.62,phase_rel_dtemp=180.0,rel_dg=10.0,phase_rel_dg=34.0):
        self.l=l
        self.m=m
        self.rel_dr=rel_dr
        self.k=k
        self.frequency=frequency
        self.phase_offset=phase_offset
        self.rel_dtemp=rel_dtemp
        self.phase_rel_dtemp=phase_rel_dtemp
        self.rel_dg=rel_dg
        self.phase_rel_dg=phase_rel_dg

#define star data
class star_data(BaseModel):
    mass:float
    radius:float
    effective_temperature:float
    vsini:float
    inclination_angle:float
    
#define phases of pulsation
class UniformTime(BaseModel):
    start:float
    end:float
    step:float

class ExplicitTime(BaseModel):
    collection:List[float]

class TimePoints(BaseModel):
    Uniform:Optional[UniformTime]=None
    Explicit:Optional[ExplicitTime]=None

#define type of meshing
class SphericalStar(BaseModel):
    theta:float
    phi:float

class Mesh(BaseModel):
    Sphere:Optional[SphericalStar]=None



class PulstarConfig(BaseModel):
    mode_data:List[mode]
    star_data:star_data
    time_points:TimePoints
    mesh:Mesh
   