
Once a filter has been created and defined, the shiftStages() methods allow the user to move
the filter frequency response up and down in QUADRATURE. This means that if you try to filter
a real signal with a quadrature bandpass, then the filter may not filter out both pos/neg
frequency components. To get around this, use the shiftStages(), but then overwrite the "tapgains"
property in the filter object with the real part of the tap gains and multiply by 2 (to compensate
for the loss that occurs when throwing out the imaginary part):
   filtobj.tapgains = 2*real(filtobj.tapgains)



