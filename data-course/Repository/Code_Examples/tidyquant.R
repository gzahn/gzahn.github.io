library(tidyverse)
library(lubridate)
library(tidyquant)
library(patchwork)
library(easystats)
theme_set(theme_minimal())

df <- tidyquant::tq_get("TSLA")

df <- df %>% 
  mutate(daily_gain=open-close) %>% 
  mutate(gain_or_loss = case_when(daily_gain <=0 ~ "Loss",
                                  daily_gain >0 ~ "Gain")) %>% 
  mutate(volatility=abs(daily_gain))

p1 <- df %>% 
  ggplot(aes(x=date,y=adjusted)) +
  geom_point() +
  labs(title = "Adjusted price")

p2 <- df %>% 
  ggplot(aes(x=date,y=close-open,color=gain_or_loss)) +
  geom_point() +
  labs(title = "Daily gains")

p1 / p2

p3 <- df %>% 
  ggplot(aes(x=date,y=volatility)) +
  geom_point(alpha=.25) + geom_smooth() +
  labs(title = "Volatility")

p1/p2/p3

p4 <- df %>% 
  ggplot(aes(x=adjusted,y=volatility)) +
  geom_point(alpha=.25) + geom_smooth() +
  labs(title = "Price vs Volatility")

(p1 + p2) / (p3 + p4)



mod <- glm(data=df,
           formula= volatility ~ (date + volume) * gain_or_loss)

step <- MASS::stepAIC(mod)$formula
mod2 <- glm(data=df,
            formula = step)

aov(mod)
summary(mod2)

modelr::add_predictions(df,mod2) %>% 
  ggplot(aes(x=date)) +
  geom_point(aes(y=volatility,color=gain_or_loss),alpha=.1) +
  geom_point(aes(y=pred),color="Black")

performance(mod2)
performance::performance_accuracy(mod2)
plot(mod2)

