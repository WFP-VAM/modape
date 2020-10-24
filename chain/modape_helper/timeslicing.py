from abc import ABC, abstractmethod
import datetime
from dateutil.relativedelta import relativedelta

class TimeSlice(ABC):

    def __init__(self):
        super().__init__()

    @property
    @abstractmethod
    def Year(self):
        pass

    @property
    @abstractmethod
    def Seqno(self):
        return self.__seqno

    @property
    @abstractmethod
    def Month(self):
        pass

    @property
    @abstractmethod
    def SliceInMonth(self):
        pass

    @abstractmethod
    def Equals(self, otherSlice):
        pass

    @abstractmethod
    def greaterThan(self, otherSlice):
        pass

    @abstractmethod
    def subtract(self, deltaSlices):
        pass

    @abstractmethod
    def add(self, deltaSlices):
        pass

    @abstractmethod
    def next(self):
        pass

    @abstractmethod
    def prev(self):
        pass

    @abstractmethod
    def getDateTimeStart(self):
        pass

    @abstractmethod
    def getDateTimeMid(self):
        pass

    @abstractmethod
    def getDateTimeEnd(self):
        pass

    @abstractmethod
    def startsBeforeDate(self, someDate):
        pass

    @abstractmethod
    def endsAfterDate(self, someDate):
        pass

    @abstractmethod
    def getFormattedDate(self, fmt):
        pass

    @abstractmethod
    def __str__(self):
        pass


class ModisInterleavedOctad(TimeSlice):

    def __init__(self, year=None, seqno=None):
        super().__init__()
        self.__year = None
        self.__seqno = None
        if isinstance(year, datetime.date):
            self.__year = year.year
            self.__seqno = (year.timetuple().tm_yday // 8) + 1
            if seqno and self.getDateTimeMid() > year:
                # when mid is after the provided date, go back one time step:
                raise NotImplementedError()
        else:
            self.__year = year
            self.__seqno = seqno
            
    @property
    def Year(self):
        return self.__year

    @property
    def Seqno(self):
        return self.__seqno

    @property
    def Month(self):
        raise NotImplementedError()

    @property
    def SliceInMonth(self):
        raise NotImplementedError()

    def Equals(self, otherSlice):
        return self.Year == otherSlice.Year and self.Seqno == otherSlice.Seqno

    def greaterThan(self, otherSlice):
        raise NotImplementedError()

    def subtract(self, deltaSlices):
        return self.add(-1 * deltaSlices)

    def add(self, deltaSlices):
        octads = (self.Year * 46) + self.Seqno + deltaSlices
        if (octads % 46) == 0:
            return ModisInterleavedOctad((octads//46)-1, 46)
        else:
            return ModisInterleavedOctad(octads//46, (octads % 46))

    def next(self):
        return self.add(1)

    def prev(self):
        return self.add(-1)

    def getDateTimeStart(self):
        return datetime.datetime.strptime('{}{}'.format(self.__year, ((self.__seqno-1) * 8) + 1), '%Y%j')

    def getDateTimeMid(self):
        raise NotImplementedError()

    def getDateTimeEnd(self):
        return self.next().getDateTimeStart() - relativedelta(seconds=1)
        
    def startsBeforeDate(self, someDate):
        return self.getDateTimeStart() < datetime.datetime(someDate.year, someDate.month, someDate.day)

    def endsAfterDate(self, someDate):
        return self.getDateTimeEnd() > datetime.datetime(someDate.year, someDate.month, someDate.day, 23, 59)

    def getFormattedDate(self, fmt):
        return self.getDateTimeStart().strftime(fmt)

    def __str__(self):
        return '{}{:03d}'.format(self.Year, ((self.Seqno - 1) * 8) + 1)
        

class Dekad(TimeSlice):

    def __init__(self, year=None, seqno=None):
        super().__init__()
        self.__year = None
        self.__seqno = None
        if isinstance(year, datetime.date):
            self.__year = year.year
            self.__seqno = ((year.month - 1) * 3) + (min(2, (year.day - 1) // 10) + 1)
            if seqno and self.getDateTimeMid() > year:
                # when mid is after the provided date, go back one time step:
                dekads = (self.Year * 36) + self.Seqno - 1;
                if (dekads % 36) == 0:
                    self.__year = (dekads // 36) - 1
                    self.__seqno = 36
                else:
                    self.__year = dekads // 36
                    self.__seqno = dekads % 36
        else:
            self.__year = year
            self.__seqno = seqno

    @property
    def Year(self):
        return self.__year

    @property
    def Seqno(self):
        return self.__seqno

    @property
    def Month(self):
        return ((self.__seqno-1)//3) + 1

    @property
    def SliceInMonth(self):
        return ((self.__seqno-1) % 3)+1;

    def deltaDekads(self, otherSlice):
        return ((self.Year - otherSlice.Year) * 36) + (self.Seqno - otherSlice.Seqno)

    def subtract(self, deltaSlices):
        return self.add(-1 * deltaSlices)

    def add(self, deltaSlices):
        dekads = (self.Year * 36) + self.Seqno + deltaSlices;
        if (dekads % 36) == 0:
            return Dekad((dekads//36)-1, 36)
        else:
            return Dekad(dekads//36, (dekads % 36))

    def next(self):
        return self.add(1)

    def prev(self):
        return self.add(-1)

    def greaterThan(self, otherSlice):
        return ((self.Year * 36) + self.Seqno) > ((otherSlice.Year * 36) + otherSlice.Seqno)

    def Equals(self, otherSlice):
        return self.Year == otherSlice.Year and self.Seqno == otherSlice.Seqno

    def getDateTimeMid(self):
        return self.getDateTimeStart() + relativedelta(days=4)

    def getDateTimeStart(self):
        return datetime.datetime(self.Year, self.Month, ((self.SliceInMonth - 1) * 10) + 1)

    def getDateTimeEnd(self):
        return self.next().getDateTimeStart() - relativedelta(seconds=1)

    def startsBeforeDate(self, someDate):
        return self.getDateTimeStart() < datetime.datetime(someDate.year, someDate.month, someDate.day)

    def endsAfterDate(self, someDate):
        return self.getDateTimeEnd() > datetime.datetime(someDate.year, someDate.month, someDate.day, 23, 59)

    def getFormattedDate(self, fmt):
        return self.getDateTimeMid().strftime(fmt)

    def __str__(self):
        return '{}{:02d}{}'.format(self.Year, self.Month, self.SliceInMonth)