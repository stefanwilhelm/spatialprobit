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
    \item{\code{adoption.date}}{an ordered factor with levels 
      \code{November, 1953} < \code{December, 1953} < \code{January, 1954} < 
      \code{February, 1954} < \code{March, 1954} < \code{April, 1954} < 
      \code{May, 1954} < \code{June, 1954} < \code{July, 1954} < 
      \code{August, 1954} < \code{September, 1954} < \code{October, 1954} < 
      \code{November, 1954} < \code{December, 1954} < 
      \code{December/January, 1954/1955} < \code{January/February, 1955} < 
      \code{February, 1955} < \code{no prescriptions found} < 
      \code{no prescription data obtained}}
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
The description of the data set from \url{http://moreno.ss.uci.edu/data.html\#ckm}:

This data set was prepared by Ron Burt. 
He dug out the 1966 data collected by Coleman, Katz and Menzel on medical innovation. 
They had collected data from physicians in four towns in Illinois, Peoria, 
Bloomington, Quincy and Galesburg. 

They were concerned with the impact of network ties on the physicians' 
adoprion of a new drug, tetracycline. 
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

The codes are:\cr
City (: 1 Peoria, 2 Bloomington, 3 Quincy, 4 Galesburg

Adoption Date:
\tabular{rl}{
1 \tab November, 1953 \cr
2 \tab December, 1953 \cr 
3 \tab January, 1954 \cr 
4 \tab February, 1954 \cr 
5 \tab March, 1954 \cr 
6 \tab April, 1954 \cr 
7 \tab May, 1954 \cr 
8 \tab June, 1954 \cr 
9 \tab July, 1954 \cr 
10 \tab August, 1954 \cr 
11 \tab September, 1954 \cr 
12 \tab October, 1954 \cr 
13 \tab November, 1954 \cr 
14 \tab December, 1954 \cr 
15 \tab December/January, 1954/1955 \cr 
16 \tab January/February, 1955 \cr 
17 \tab February, 1955 \cr 
18 \tab no prescriptions found \cr 
98 \tab no prescription data obtained 
}

Year started in the profession
\tabular{rl}{
1 \tab 1919 or before \cr 
2 \tab 1920-1929 \cr 
3 \tab 1930-1934 \cr 
4 \tab 1935-1939 \cr 
5 \tab 1940-1944 \cr 
6 \tab 1945 or later \cr 
9 \tab no answer 
}

Have you attended any national, regional or state conventions of professional 
societies during the last 12 months? [if yes] Which ones? 
\tabular{rl}{
0 \tab none \cr 
1 \tab only general meetings \cr 
2 \tab specialty meetings \cr 
9 \tab no answer \cr
}

Which medical journals do you receive regularly?
\tabular{rl}{ 
1 \tab two \cr 
2 \tab three \cr 
3 \tab four \cr 
4 \tab five \cr 
5 \tab six \cr 
6 \tab seven \cr 
7 \tab eight \cr 
8 \tab nine or more \cr
9 \tab no answer 
}

With whom do you actually spend more of your free time -- doctors or non-doctors? 
\tabular{rl}{
1 \tab non-doctors \cr 
2 \tab about evenly split between them \cr 
3 \tab doctors \cr 
9 \tab mssing; no answer, don't know 
}

When you are with other doctors socially, do you like to talk about medical matter? 
\tabular{rl}{
1 \tab no \cr 
2 \tab yes \cr 
3 \tab don't care \cr 
9 \tab missing; no answer, don't know 
}

Do you belong to any club or hobby composed mostly of doctors?
\tabular{rl}{ 
0 \tab no \cr 
1 \tab yes \cr 
9 \tab no answer \cr
} 

Would you tell me who are your three friends whom you see most often socially? What is [their] occupation?
\tabular{rl}{
1 \tab none are doctors \cr 
2 \tab one is a doctor \cr 
3 \tab two are doctors \cr 
4 \tab three are doctors \cr 
9 \tab no answer 
}

How long have you been practicing in this community? 
\tabular{rl}{
1 \tab a year or less \cr 
2 \tab more than a year, up to two years \cr 
3 \tab more than two years, up to five years \cr 
4 \tab more than five years, up to ten years \cr 
5 \tab more than ten years, up to twenty years \cr 
6 \tab more than twenty years \cr 
9 \tab no answer 
}

About how many office visits would you say you have during the average week at this time of year?
\tabular{rl}{
1 \tab 25 or less \cr 
2 \tab 26-50 \cr 
3 \tab 51-75 \cr 
4 \tab 76-100 \cr 
5 \tab 101-150 \cr 
6 \tab 151 or more \cr
9 \tab missing; no answer, don't know
} 

Are there other physicians in this building? 
[if yes] Other physicians in same office or with same waiting room? 
\tabular{rl}{
1 \tab none in building \cr 
2 \tab some in building, but none share his office or waiting room \cr 
3 \tab some in building sharing his office or waiting room \cr 
4 \tab some in building perhaps sharing his office or waiting room \cr 
9 \tab no answer 
}

Do you specialize in any particular field of medicine? [if yes] What is it?
\tabular{rl}{
1 \tab GP, general practitioner \cr 
2 \tab internist \cr 
3 \tab pediatrician \cr 
4 \tab other specialty \cr 
9 \tab no answer
}
}

\source{
The data set is reproduced from \url{http://moreno.ss.uci.edu/data.html\#ckm}
with the friendly permission of Prof. Lin Freeman.
}

\references{
Burt, R. (1987). Social contagion and innovation: Cohesion versus structural equivalence. \emph{American Journal of Sociology}, \bold{92}, 1287--1335. 

Coleman, James, Elihu Katz and Herbert Menzel (1957). The Diffusion of an Innovation Among Physicians, \emph{Sociometry}, \bold{20}, 253--270. 

Coleman, J.S., E. Katz, and H. Menzel (1966). \emph{Medical Innovation: A Diffusion Study}. New York: Bobbs Merrill. 

Valente, T. W. (1995). \emph{Network Models of the Diffusion of Innovations}. Cresskill, NJ: Hampton Press. 

Van den Bulte, C. and G. L. Lilien. (2001). Medical Innovation Revisited: Social Contagion versus Marketing Effort, 
  \emph{American Journal of Sociology}, \bold{106}, 1409--1435. 
}
\examples{
data(CKM)
}
\keyword{datasets}
