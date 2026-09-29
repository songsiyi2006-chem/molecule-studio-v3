"""V3 catalog routes, shared by the offline browser and external clients."""
from typing import Annotated, Literal
import sqlite3
from fastapi import APIRouter, HTTPException, Query
from .catalog import Catalog, QuantumCatalog
router=APIRouter(tags=['Reference catalogs'])
def checked(call):
    try:return call()
    except ValueError as e:raise HTTPException(422,str(e)) from e
    except (FileNotFoundError, OSError, sqlite3.Error) as e:raise HTTPException(503,'参考数据库或索引暂不可用，请运行数据重建脚本后重试。') from e

@router.get('/database/advanced')
def advanced(q:Annotated[str,Query(max_length=2000)]='',mode:Literal['text','exact','similarity']='text',source:Annotated[str|None,Query(max_length=50)]=None,target:Literal['logS','logD74','hydration_free_energy']|None=None,limit:Annotated[int,Query(ge=1,le=100)]=20,offset:Annotated[int,Query(ge=0,le=1000000)]=0,threshold:Annotated[float,Query(ge=0,le=1)]=.2):
    return checked(lambda:Catalog().search(q,mode,source,target,limit,offset,threshold))

@router.get('/database/detail/{molecule_id}')
def detail(molecule_id:int):
    result=checked(lambda:Catalog().detail(molecule_id))
    if result is None:raise HTTPException(404,'未找到该分子。')
    return result

@router.get('/quantum/stats')
def quantum_stats():return checked(lambda:QuantumCatalog().stats())

@router.get('/quantum/search')
def quantum_search(q:Annotated[str,Query(max_length=2000)]='',limit:Annotated[int,Query(ge=1,le=100)]=20,offset:Annotated[int,Query(ge=0,le=1000000)]=0):
    return checked(lambda:QuantumCatalog().search(q,limit,offset))
