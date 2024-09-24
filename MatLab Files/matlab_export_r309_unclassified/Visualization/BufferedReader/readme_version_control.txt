Jeff H.  11-13-2017
 - accumulated versions of bufferedReader.m. This script is the LATEST.
 -  Some of the updates:
    - made jumpToPage() FAST by modifying shiftPage() to take a 3rd argument that is a switch for either loading data
      or only updating the file pointers.
      ** included constant class properties "enum_fast" and "enum_load" flags that are used in the updated shiftPage()
    - partialShiftPage is still broken, but has been included in this latest version for future work
    - added property "headBlockBytes" to class to account for the initial file header
