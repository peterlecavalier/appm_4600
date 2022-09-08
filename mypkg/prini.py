class prini:
    def __init__(self,xtype,text,x):


      self.type = xtype
      self.text = text
      self.val = x
      
    def print(self):
    
      if self.type == 'real':
        print(self.text, '%16.16e' % self.val)
      elif self.type == 'inter':
         print(self.text, '%d' % self.val)  
      elif self.type == 'comp':
         print(self.text, '%0.4f + %0.4fi' % (self.val.real, self.val.imag))
      else:
        print('type not supported')
        
