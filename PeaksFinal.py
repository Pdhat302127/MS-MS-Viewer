   
from base64 import b64decode
from array import array
import xml.etree.ElementTree as ET
import sys

def peak_extraction(afile, ns, scan):
    name_space=False
    scan_number=False
    try:
        xml_document = ET.iterparse(afile)
        for event, element in xml_document:
            if element.tag == ns + 'scan':
                name_space=True
                if element.attrib.get('num') == scan:
                    scan_number=True
                    peak = element.find(ns + 'peaks')
                    
                    if peak is not None:
                        peaks = array('f', b64decode(peak.text))
                        if sys.byteorder != 'big':
                            peaks.byteswap()
                
                        mzs = peaks[::2]  
                        ints = peaks[1::2]  
                        return mzs, ints  
                element.clear()
                
        if name_space is not True:
            raise ValueError("ns:"+ ns + " " + "might have an error or not defined")        
        if scan_number is not True:
            raise ValueError("Scan:"+ scan + " " + "not found in xml file")
        

    except Exception as e:
        raise ValueError("An error occurred:" + " " + str(e))
    
      
