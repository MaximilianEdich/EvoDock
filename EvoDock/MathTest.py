import matplotlib.pyplot as plt

S = 0.9
I = 1 - S
R = 0

a = 0.25
µ = 0.125

runTime = 200

Ss = []
Is = []
Rs = []

for i in range(runTime):
    newS = S - a * S * I
    newI = I + a * S * I - µ * I
    newR = R + µ * I
    S = newS
    I = newI
    R = newR
    Ss.append(S)
    Is.append(I)
    Rs.append(R)

print(S)
print(I)
print(R)
print(µ / a)

# Plot results
plt.plot(Ss)
plt.plot(Is)
plt.plot(Rs)
plt.show()