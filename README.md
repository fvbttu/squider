# squider
An ODE model for COVID-19 infection and death rates

SQUIDER is an SIR-based compartmental model designed to take into account certain features relevant to COVID-19:
 - power law incidence rate to more accurately depict compartment interaction in jurisdictions with heterogeneous population density,
 - separate compartments for detected and undetected infected and recovered populations (silent spreaders),
 - time-dependent applications and relaxations of restrictions on social contact (e.g. stay-at-home orders),
 - quarantine / hospital isolation of confirmed infected individuals,
 - possible loss of immunity for recovered individuals in the short-to-medium time scale (as is the case with other coronaviruses, such as the common cold).

### usage
The model has been implemented in Matlab (see code). To obtain values for parameters fits are done to case-rate and death time series such as are available from the [Johns Hopkins University's repository](https://github.com/CSSEGISandData/COVID-19) on GitHub. This may be preceeded by smoothing spline fits to estimate the number of separate interventions that may have occurred. Once a good fit has been made, the parameters can be used to obtain projections for future dates by running the ODE solver on those parameters to the day in question, with any hypothetical interventions or relaxations beyond the last day of fit data appended to current parameter set.

### people
SQUIDER comes out of the **Hussain Lab at Texas Tech University** COVID-19 working group:
- Fazle Hussain, Ph.D., President's Endowed Distinguished Chair in Engineering, Science, & Medicine, TTU
  * (email: `fazle.hussain@ttu.edu`)
- Zeina Khan, Ph.D., Research Assistant Professor, TTU Mechanical Engineering (email: `zeina.khan@ttu.edu`)
- Frank Van Bussel,  Ph.D. (email: `frank.van-bussel@ttu.edu`)

### preprint
A [write-up of our current results](https://arxiv.org/abs/2006.05955) is available on Arxiv.
