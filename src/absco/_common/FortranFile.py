# for Python 3 compatibility
from __future__ import print_function

import struct
import types

class FortranFile:
    sizeFormat='I'

    def __init__(self,fileName,write=False,network=False):
        if write:
            self.file=open(fileName,'wb')
        else:
            self.file=open(fileName,'rb')

        self.write=write
    
        if network: self.byteOrder='!'
        else: self.byteOrder='='
            
    def writeRecord(self,vector):
        pass
    
    def readFloatVector(self):
        return self.readRecord('f')

    def readDoubleVector(self):
        return self.readRecord('d')

    def readFormatVector(self,formatString):
        return self.readRecord(formatString)

    def readFormatData(self,formatString):
        record=self.getRecord()
        if record:
            vectorFmt='%s%s'%(self.byteOrder,formatString)
            try:
                return struct.unpack(vectorFmt,record)
            except:
                return None
        return None

    def getRecord(self):
        stringSize = struct.calcsize(self.sizeFormat)
        try:
            size = struct.unpack(self.sizeFormat, self.file.read(stringSize))[0]
            buff = self.file.read(size)
            otherSize = struct.unpack(self.sizeFormat, self.file.read(stringSize))[0]
            if otherSize != size: return None
            return buff
        except:
            return None
        
    def readRecord(self,id):
        record=self.getRecord()
        if record:
            size=len(record)/struct.calcsize(id)
            vectorFmt='%s%0d%s'%(self.byteOrder,size,id)
            try:
                return struct.unpack(vectorFmt,record)
            except:
                return None
        return None
    
    def reset(self):
        self.file.seek(0)

    def writeSpecialFormatVector(self,data,formatString,outputSize):
        buffer=''
        for index,value in enumerate(data):
            print(formatString[index],value)
            if type(value) in [types.ListType,types.TupleType]:
                for x in value: buffer += struct.pack(formatString[index],x)
            else:
                buffer += struct.pack(formatString[index],value)
        for index in range(len(buffer),outputSize):buffer+=struct.pack('b',0)

        size = struct.pack(self.sizeFormat,(len(buffer)))
        print(len(buffer))
        self.file.write(size+buffer+size)

    def writeFormatVector(self,data,formatString,outputSize=None):
        buffer=''
        for index,value in enumerate(data):
            if type(value) in [types.ListType,types.TupleType]:
                for x in value: buffer += struct.pack(formatString[index],x)
            else:
                buffer += struct.pack(formatString[index],value)

        if outputSize is not None:
            for index in range(len(buffer),outputSize):buffer+=struct.pack('b',0)

        size = struct.pack(self.sizeFormat,(len(buffer)))
        self.file.write(size+buffer+size)

    def writeFloatVector(self,data):
        formatString=''
        for i in data: formatString+='f'
        self.writeFormatVector(data,formatString)
        
    def close(self):
        if self.file:
            self.file.close()
            self.file=None

    def __del__(self):
        self.close()

def readReflectance(fileName):
    fortranFile=FortranFile(fileName)
    format='114d2i5d4i10d'
    data=fortranFile.readFormatData(format)
    #print data

    OK=True
    length=0
    wn1=None
    refls=[]
    while OK:
        try:
            v1,v2,dv,nWN=fortranFile.readFormatVector('ddfi')
            #print v1,v2,dv,nWN
        except:
            nWN=-1
            
        if nWN>0:
            data=fortranFile.readFloatVector()
            refls+=data
            if wn1 is None:
                wn1=v1
                wndv=dv
            wn2=v2
        else: OK=False

    return wn1,wn2,wndv,refls

def writeReflectance(fileName,v1,v2,dv,reflCoeff):
    def refl(wn,coeff):
        return coeff[0]+coeff[1]*wn+coeff[2]*wn**2
    
    maxWNs=2400
    nFillHeader=264

    fillHeader=map(lambda x:0.,range(nFillHeader))

    fortranFile=FortranFile(fileName,write=True)
    fortranFile.writeFloatVector(fillHeader)

    nWNs=int(round((v2-v1)/dv)+1)
    if nWNs>maxWNs:
        dv=(wn2-wn1)/(maxWNs-1)
        nWNs=int(round((v2-v1)/dv)+1)
    
    fortranFile.writeFormatVector((v1,v2,dv,nWNs),'ddfi')
    refl=map(lambda x:refl(v1+x*dv,reflCoeff),range(nWNs))
    fortranFile.writeFloatVector(refl)
    fortranFile.writeFormatVector((v2,v2,dv,-99),'ddfi')
    fortranFile.close()

if __name__ == '__main__':
    fileName='SOL.REFLECTANCE'
    wn1,wn2,wndv,data=readReflectance(fileName)
    print(wn1,wn2,wndv,len(data),min(data),max(data))

