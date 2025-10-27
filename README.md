# Met4H2 software framework for uncertainty calculations
This repository contains a software platform implementing a methodology for the 
calculation of the total quantity and total energy along the hydrogen supply chain. 
The software allows the specification of input parameters and associated uncertainties,
selection of measurands, calculation of combined uncertainties, and reporting 
functionality of uncertainty budgets. 

The software has been developed in a modular fashion such that it can be expanded 
with, e.g., new calculation models for auxiliary parameters, novel measurement methods 
and equipment, etc.

The uncertainty evaluation is based on ISO/IEC Guide 98 (Guide to the expression of Uncertainty 
in Measurement) and it has been adapted for specific use regarding combined measurmeent 
uncertainty of quantities and energy content of energy gases such as hydrogen-based fluids.

## Dependencies Requirements
The scripts have been developed using Python 3.9.17.

Please refer to the *venvm4h2.yml* file for more details on the required libraries.

## Documentation
Please refer to the *m4h2_doc.pdf* file for the software documentation.

## Related projects
EPM [21GRD05 Met4H2](https://met4h2.eu/) - Metrology for the hydrogen supply chain (Oct 2022 - Sep 2025).

## Acknowledgment
The project 21GRD05 (Met4H2) has received funding from the European Partnership on Metrology, 
co-financed from the European Union's Horizon Europe Research and Innovation Programme and by 
the Participating States.

## License
This code is provided under the Creative Commons Attribution 4.0 International Public License.
