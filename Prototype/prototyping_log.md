
# 1st Version - [✅]
> System
- event
  - types: interested / not interested
  - each event single value
  - distribution - normal
- triggering system
  - above threshold - interested and capture
  - below threshold - not interested and ignore

> Estimator
- assume we know the ignored ones
- use ground truth for feedback

> Controller & Planner
- capture an uninterested one, threshold goes up a little bit
- capture an interested one, threshold goes down a little bit

![illustration](Prototype_V1/SATM_Prototype_V1.png)

> Modify the Parameters to Control

- the speed to converge
- the fluctuation
- etc

# 2nd Version - [✅]
- Feature - Metrics Calculation on the fly

> Results

![illustration](Prototype_V2/SATM_Prototype_V2.png)

> The control strategy is the same as the 1st version
> Online metrics calculation is implemented

![illustration](Prototype_V2/SATM_Prototype_Metrics_V2.png)

# 3rd Version - [under development]
- New Feature - Metrics Guided Control

![illustration](Prototype_V3/SATM_Prototype_V3.png)
- stage 1: yes or no feedback
- stage 2: metrics guided control - not working well now

# Ideas：
- traditional controller
- magnitude & duration

the idea is similar to the pinoeer sensing to get knowledge to properly setup the thresholds for magnitude and duration
# environment-aware control framework - know the environment first then control; to tackle with the unseen events, plus human knowledge

Part 1:
Sense All and Clustering - threshold
- extant classes
- new classes

Part 2:
Judge whether the event is interested or not  - (1) criterions (2) known classes

Part 3:
based on the context 
