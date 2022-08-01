# Function to find pseudo values. (Sözde Değerlerin Fonksiyonu)
# Jackknife yöntemi ile elde ediliyor.

# NOT: Sözde değerler hesaplanırken tek reader ve test üzerinden hesaplama yapılmaktadır.
#      Veriler MRMC modunda ise sözde değerler hesaplanmadan önce verinin bütün test:reader kombinasyonları için split edilmesi
#      ve analizlerin her alt grup için ayrı ayrı yapılması gerekmektedir.
pseudo_values <- function(object, statusVariable, marker, event = "1", higherValuesDiseased = TRUE){
  # Args:
  #   object: data.frame veya matris. Analiz edilecek olan veri setidir.
  #   statusVariable: status değişkeninin adı girilmelidir.
  #   marker: 
  status <- object[ ,statusVariable]
  marker <- object[ ,marker]
  
  rocRes <- rocdata(status = status, marker = marker, 
                    event = event, higherValuesDiseased = higherValuesDiseased)
  
  AUC <- as.numeric(rocRes)
  
  pseudo <- NULL
  for (i in 1:length(marker)){
    status.i <- status[-i]
    marker.i <- marker[-i]
    
    rocRes.i <- rocdata(status = status.i, marker = marker.i, 
                        event = event, higherValuesDiseased = higherValuesDiseased)
    
    AUC.i <- as.numeric(rocRes.i)
    ## Pseudo değerin jacknife değeri.
    pseudo.i <- length(status)*AUC - (length(status) - 1)*AUC.i
    pseudo <- c(pseudo, pseudo.i)
  }
  
  object <- data.frame(object, pseudo = round(pseudo, 8))
  return(object)
}