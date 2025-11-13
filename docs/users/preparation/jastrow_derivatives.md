---
layout: default
title: Jastrow Derivatives
nav_order: 11
parent: Input files
authors:
    - Ravindra Shinde
tags:
    - CHAMP
    - jastrow derivative
---

# Jastrow derivatives file
The Jastrow derivative parameters can be provided using this file. It has the following format [Example: water].

```python
jasderiv
4 4 5 15 15 0 0 nparma,nparmb,nparmc,nparmf
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
  3 4 5 6 (iwjasa(iparm),iparm=1,nparma)
2 3 4 5 6 (iwjasb(iparm),iparm=1,nparmb)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
3 5 7 8 9         11 13 14 15 16     17 18 20 21 23 (c(iparmj),iparmj=1,nparmc)
end
```