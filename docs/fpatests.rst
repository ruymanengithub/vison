.. _fpa-tests:

FPA Tests
=========

.. currentmodule:: vison.fpatests

These are tests made at FPA level, instead of block-level. The core of this is the data models in:

#. :ref:`fpa_dm.py`
#. vison/fpa/fpa.py

There are only a few tests defined (they are defined in :ref:`CEAFPAcampaign`), all for **single images**:

#. FWD_WARM: forward acquisition in warm. Analysis of dark current ramps.
#. CHINJ: basic test of charge injection (in cold).
#. DARK: basic test of a long exposure dark, in cold.
#. BIAS_RWDVS_WARM: for bias frames acquired in warm, with back clocking in vertical and serial (RWDVS).
#. BIAS_RWDV_WARM: for bias frames acquired in warm, with back clocking only in vertical direction (RWDV).
#. BIAS_RWDVS_COLD: same as BIAS_RWDVS_WARM but for detectors in cold.
#. BIAS_RWDV_WARM: same as BIAS_RWDV_WARM but for detectors in cold.
#. BIAS_FWD_COLD: for bias frames acquired in cold, forward clocking.




fpamaster.py
------------

.. automodule:: vison.fpatests.fpamaster
        :members:


fpatask
-------

.. automodule:: vison.fpatests.fpatask
        :members:
