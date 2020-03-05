Class17: Corona virus investigation
================
Aries Chavira
3/4/2020

Here we analyze infection data for the 2019 novel Coronavirus COVID-19 (2019-nCoV) epidemic. The raw data is pulled from the Johns Hopkins University Center for Systems Science and Engineering (JHU CCSE) Coronavirus repository.

A CSV file is available here <https://github.com/RamiKrispin/coronavirus-csv>

``` r
url <- "https://tinyurl.com/COVID-2019"
virus <- read.csv(url)

tail(virus)
```

    ##      Province.State Country.Region     Lat     Long       date cases      type
    ## 2675         Shanxi Mainland China 37.5777 112.2922 2020-03-03     5 recovered
    ## 2676        Sichuan Mainland China 30.6171 102.7103 2020-03-03     8 recovered
    ## 2677        Tianjin Mainland China 39.3054 117.3230 2020-03-03    13 recovered
    ## 2678       Xinjiang Mainland China 41.1129  85.2401 2020-03-03     2 recovered
    ## 2679         Yunnan Mainland China 24.9740 101.4870 2020-03-03     1 recovered
    ## 2680       Zhejiang Mainland China 29.1832 120.0934 2020-03-03    24 recovered

Questions to answer 1. How many total infected cases are there 2. How many deaths linked to infection have there been 3. What is the overall death rate of the disease of the entire dataset 4. What is the death rate in China respectively 5. What is the death rate in Italy, US, and Iran?

1.  How many total infected cases are there

``` r
total_cases <- sum(virus$cases)
```

1.  How many deaths linked to infection have there been

``` r
death <- virus$type=="death"
death_cases <- sum(virus[death,"cases"])
```

1.  What is the overall death rate of the disease of the entire dataset

``` r
3160/144233 *100
```

    ## [1] 2.190899

1.  What is the death rate in China respectively

``` r
main_china <- virus$Country.Region=="Mainland China"
main_china <- virus[main_china, ]

china_dead <- main_china$type=="death"
dr_china <- sum(main_china[china_dead, "cases"])
dr_china 
```

    ## [1] 2945

``` r
total_china <- sum(main_china$cases)
dr_in_china <- dr_china / total_china * 100
dr_in_china
```

    ## [1] 2.256705

> Death rate in ITALY

``` r
italy <- virus$Country.Region=="Italy"
total_italy <- virus[italy, ]

italy_deaths <- total_italy$type=="death"
dr_italy <- sum(total_italy[italy_deaths, "cases"])
dr_italy
```

    ## [1] 79

``` r
total_infected_italy <- sum(total_italy$cases)
dr_in_italy <- dr_italy/ total_infected_italy * 100
dr_in_italy
```

    ## [1] 2.88216

> Death rate in IRAN

``` r
iran <- virus$Country.Region=="Iran"
only_iran <- virus[iran, ]

iran_dead <- only_iran$type=="death"
dr_iran <- sum(only_iran[iran_dead, "cases"])
dr_iran
```

    ## [1] 77

``` r
total_infected_iran <- sum(only_iran$cases)
dr_in_iran <- (dr_iran / total_infected_iran * 100)
dr_in_iran
```

    ## [1] 2.847633

> Death rate in the US

``` r
us <- virus$Country.Region=="US"
total_us <- virus[us, ]

us_deaths <- total_us$type=="death"
dr_us <- sum(total_us[us_deaths, "cases"])
dr_us
```

    ## [1] 7

``` r
cases_us <- sum(total_us$cases)
cases_us
```

    ## [1] 137

``` r
dr_in_us <- (dr_us / cases_us * 100)
dr_in_us
```

    ## [1] 5.109489

``` r
dr_china
```

    ## [1] 2945

``` r
dr_italy
```

    ## [1] 79

``` r
dr_iran
```

    ## [1] 77

``` r
dr_in_china
```

    ## [1] 2.256705

``` r
dr_in_iran
```

    ## [1] 2.847633

``` r
dr_in_italy
```

    ## [1] 2.88216

``` r
dr_in_us
```

    ## [1] 5.109489
