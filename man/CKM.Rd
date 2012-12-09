\name{CKM}
\alias{CKM}
\alias{A1}
\alias{A2}
\alias{A3}
\alias{W1}
\alias{W2}
\alias{W3}
\docType{data}
\title{Coleman, Katz, Menzel "Innovation among Physicians" dataset}
\description{
 The classic Coleman's Drug Adoption dataset "Innovation among Physicians"
 for studying the information diffusion through social networks.
}
\usage{data(CKM)}
\format{
  A data frame \code{CKM} with 246 observations on the following 13 variables.
  \describe{
    \item{\code{city}}{a numeric vector; City: 1 Peoria, 2 Bloomington, 3 Quincy, 4 Galesburg}
    \item{\code{adoption.date}}{an ordered factor with levels \code{November, 1953} < \code{December, 1953} < \code{January, 1954} < \code{February, 1954} < \code{March, 1954} < \code{April, 1954} < \code{May, 1954} < \code{June, 1954} < \code{July, 1954} < \code{August, 1954} < \code{September, 1954} < \code{October, 1954} < \code{November, 1954} < \code{December, 1954} < \code{December/January, 1954/1955} < \code{January/February, 1955} < \code{February, 1955} < \code{no prescriptions found} < \code{no prescription data obtained} < \code{no prescription data obtained}}
    \item{\code{med_sch_yr}}{years in practice}
    \item{\code{meetings}}{meetings attended}
    \item{\code{jours}}{journal subscriptions}
    \item{\code{free_time}}{free time activities}
    \item{\code{discuss}}{discussions}
    \item{\code{clubs}}{club memberships}
    \item{\code{friends}}{friends}
    \item{\code{community}}{time in the community}
    \item{\code{patients}}{patient load}
    \item{\code{proximity}}{physical proximity to other physicians}
    \item{\code{specialty}}{medical specialty}
  }
  
  Three 246 \eqn{\times}{x} 246 binary peer matrices \code{A1},\code{A2},\code{A3} 
  for three different social relationships/networks: "Advice", "Discussion", "Friend".
  
  Three 246 \eqn{\times}{x} 246 spatial weight matrices \code{W1}, \code{W2} and \code{W3} 
  from built from adjacency matrices \code{A1},\code{A2},\code{A3}.
}
\details{
This data set was prepared by Ron Burt. 
He dug out the 1966 data collected by Coleman, Katz and Menzel on medical innovation. 
They had collected data from physicians in four towns in Illinois, Peoria, Bloomington, Quincy and Galesburg. 

They were concerned with the impact of network ties on the physicians' adoprion of a new drug, tetracycline. 
Three sociometric matrices were generated. 
One was based on the replies to a question, 
"When you need information or advice about questions of therapy where do you usually turn?" 
A second stemmed from the question 
"And who are the three or four physicians with whom you most often find yourself discussing cases or therapy in the course of an ordinary week -- last week for instance?" 
And the third was simply 
"Would you tell me the first names of your three friends whom you see most often socially?"

In addition, records of prescriptions were reviewed and a great many other questions were asked. 
In the CKM data I have included 13 items: city of practice, recorded date of tetracycline adoption date, 
years in practice, meetings attended, journal subscriptions, free time activities, discussions, club memberships, 
friends, time in the community, patient load, physical proximity to other physicians and medical specialty.

The codes are:
City (: 1 Peoria, 2 Bloomington, 3 Quincy, 4 Galesburg

Adoption Date:
1 November, 1953 
2 December, 1953 
3 January, 1954 
4 February, 1954 
5 March, 1954 
6 April, 1954 
7 May, 1954 
8 June, 1954 
9 July, 1954 
) 10 August, 1954 
11 September, 1954 
12 October, 1954 
13 November, 1954 
14 December, 1954 
15 December/January, 1954/1955 
16 January/February, 1955 
17 February, 1955 
18 no prescriptions found 
98 no prescription data obtained 

Year started in the profession
1 1919 or before 
2 1920-1929 
3 1930-1934 
4 1935-1939 
5 1940-1944 
6 1945 or later 
9 no answer 

Have you attended any national, regional or state conventions of professional societies during the last 12 months? [if yes] Which ones? 
0 none 
1 only general meetings 
2 specialty meetings 
9 no answer

Which medical journals do you receive regularly? 
1 two 
2 three 
3 four 
4 five 
5 six 
6 seven 
7 eight 
8 nine or more
9 no answer 

With whom do you actually spend more of your free time -- doctors or non-doctors? 
1 non-doctors 
2 about evenly split between them 
3 doctors 
9 mssing; no answer, don't know 

When you are with other doctors socially, do you like to talk about medical matter? 
1 no 
2 yes 
3 don't care 
9 missing; no answer, don't know 

Do you belong to any club or hobby composed mostly of doctors? 
0 no 
1 yes 
9 no answer 

Would you tell me who are your three friends whom you see most often socially? What is [their] occupation?
1 none are doctors 
2 one is a doctor 
3 two are doctors 
4 three are doctors 
9 no answer 

How long have you been practicing in this community? 
1 a year or less 
2 more than a year, up to two years 
3 more than two years, up to five years 
4 more than five years, up to ten years 
5 more than ten years, up to twenty years 
6 more than twenty years 
9 no answer 

About how many office visits would you say you have during the average week at this time of year?
1 25 or less 
2 26-50 
3 51-75 
4 76-100 
5 101-150 
6 151 or more
9 missing; no answer, don't know 

Are there other physicians in this building? [if yes] Other physicians in same office or with same waiting room? 
1 none in building 
2 some in building, but none share his office or waiting room 
3 some in building sharing his office or waiting room 
4 some in building perhaps sharing his office or waiting room 
9 no answer 

Do you specialize in any particular field of medicine? [if yes] What is it? 
1 GP, general practitioner 
2 internist 
3 pediatrician 
4 other specialty 
9 no answer
}

\source{
The data set is reproduced from http://moreno.ss.uci.edu/data.html\#ckm 
with the friendly permission of Prof. Lin Freeman.
}

\references{
Burt, Ronald S. 1987. "Social Contagion and Innovation: Cohesion Versus Structural Equivalence," AJS 92: 1287-1335. 

Burt, R., (1987). Social contagion and innovation: Cohesion versus structural equivalence. American Journal of Sociology, 1987. 92: p. 1287-1335. 

Coleman, James, Elihu Katz and Herbert Menzel. 1957. "The Diffusion of an Innovation Among Physicians," Sociometry, 20:253-270. 

Coleman, J.S., E. Katz, and H. Menzel, 1966. Medical Innovation: A Diffusion Study. New York: Bobbs Merrill. 

Valente, T. W. (1995). Network Models of the Diffusion of Innovations. Cresskill, NJ: Hampton Press. 

Van den Bulte, C. and G. L. Lilien. 2001. "Medical Innovation Revisited: Social Contagion versus Marketing Effort,"AJS 106: 1409-1435. 
}
\examples{
data(CKM)
}
\keyword{datasets}
